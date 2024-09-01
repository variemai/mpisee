/* (C) Jesper Larsson Traff, January 2018, May 2018, October 2018 */
/* Quicksort/HyperQuicksort without pairwise element exchange
   (compared against standard parallel Quicksort and HyperQuicksort) 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

#include <assert.h>

#include <mpi.h>

#define inline

#define PARQSORT
#define HYPQSORT
#define TUWQSORT
#define WUTQSORT
//#define ITEQSORT

#define COMSORT // combine exchange and exchange-free variants
#define PARFACT 1500
#define HYPFACT 1500

#define USEGROUPS        // Groups for communicator splitting
//#define USEGROUPCREATE   // Groupwise collective communicator creation
//#define USEBLOCKING      // Blocking communication
#define USEMEMCPY        // Blocking memcpy
//#define USESENTINEL      // Sentinels in partitioning
#define USEINTERPOLATION // Average over local pivots
//#define STANDARDMERGE    // Textbook merge

#define RANGE 500     // for random indicator selection
                      // (thanks to Francesco Versaci)
#define K 3           // global median of K (for now, K=1, or K=3)
#define MMS 10        // min-max sample size

#define SEQREPEAT 5
#define REPEAT 43
#define WARMUP 5
#define MICRO  1000000.0

// Inputtypes
#define UNI 0
#define ASC 1
#define PRM 2

//typedef int base_t;
//#define BASE_T MPI_INT
//#define BASEINT_T MPI_2INT
//#define MAXBASE INT_MAX
//typedef long base_t;
//#define BASE_T MPI_LONG
//#define BASEINT_T MPI_LONG_INT
//#define MAXBASE LONG_MAX
//typedef float base_t;
//#define BASE_T MPI_FLOAT
//#define BASEINT_T MPI_FLOAT_INT
//#define MAXBASE FLT_MAX
typedef double base_t;
#define BASE_T MPI_DOUBLE
#define BASEINT_T MPI_DOUBLE_INT
#define MAXBASE DBL_MAX

typedef struct {
  base_t val;
  int loc;
} baseint_t;

typedef struct {
  base_t val;
  base_t loc;
} basebase_t;

#define QTAG 7777

// forward declarations
base_t *tuw_qsort(base_t *a, const int n, int *m, MPI_Comm comm);
base_t *wut_merge(base_t *a, const int n, int *m, MPI_Comm comm);

// O(log w) implementation, this one is not branch-free
inline unsigned int ilog2(unsigned int n) 
{
  unsigned int d = 0;

  if ((n&0xFFFF0000)!=0x0) {
    d += 16; n >>= 16;
  }
  if ((n&0xFF00)!=0x0) {
    d += 8; n >>= 8;
  }
  if ((n&0xF0)!=0x0) {
    d += 4; n >>= 4;
  }
  if ((n&0xC)!=0x0) {
    d += 2; n >>= 2;
  }
  if ((n&0x2)!=0x0) {
    d += 1; //n >>= 1;
  }

  return d;
}

void tuw_selsum(void *invec, void *inoutvec, int *len, MPI_Datatype *type)
{
  int i;

  assert(*type==BASEINT_T);
  
  baseint_t *in    = (baseint_t*)invec;
  baseint_t *inout = (baseint_t*)inoutvec;

  for (i=0; i<*len; i++) {
    if (inout[i].loc==0) {
      inout[i].val = in[i].val;
      inout[i].loc = in[i].loc;
    } else if (in[i].loc>0) {
      inout[i].val = in[i].val+inout[i].val;
      inout[i].loc = in[i].loc+inout[i].loc;
    }
  }
}

// segmented replace: bcast
void tuw_segrep(void *invec, void *inoutvec, int *len, MPI_Datatype *type)
{
  int i;

  assert(*type==BASEINT_T);
  
  baseint_t *in    = (baseint_t*)invec;
  baseint_t *inout = (baseint_t*)inoutvec;

  for (i=0; i<*len; i++) {
    if (inout[i].loc==0) {
      inout[i].val = in[i].val;
      inout[i].loc = in[i].loc;
    }
  }
}

// select maximum indicated element
void tuw_indmax(void *invec, void *inoutvec, int *len, MPI_Datatype *type)
{
  int i;

  assert(*type==BASEINT_T);
  
  baseint_t *in    = (baseint_t*)invec;
  baseint_t *inout = (baseint_t*)inoutvec;

  for (i=0; i<*len; i++) {
    if (in[i].loc>inout[i].loc) {
      inout[i].val = in[i].val;
      inout[i].loc = in[i].loc;
    }
  }
}

MPI_Op TUW_CHOOSE;
MPI_Op TUW_SEGREP; // the segmented broadcast

#define NULLTAG 100
#define BCSTTAG 200
int bcast_segment(void *buffer, int count, MPI_Datatype datatype,
		  int isroot, MPI_Comm comm)
{
  int rank, size;
  int succ, pred;
  int d;

  void *sendbuf, *recvbuf, *nullbuf;
  int sendcount, recvcount;
  int sendtag;

  MPI_Aint lb, extent;

  MPI_Status status;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  MPI_Type_get_extent(datatype,&lb,&extent);
  nullbuf = (void*)malloc(count*extent);

  sendbuf   = buffer;
  sendcount = 0;
  sendtag   = NULLTAG;
  recvbuf   = buffer;
  recvcount = count;

  // Hillis-Steele like implementation for small problems
  for (d=1; d<size; d<<=1) {
    succ = rank+d; if (succ>=size) succ = MPI_PROC_NULL;
    pred = rank-d; if (pred<0)     pred = MPI_PROC_NULL;

    if (isroot) {
      sendcount = count;
      sendtag   = BCSTTAG;
      recvbuf   = nullbuf;
    }

    MPI_Sendrecv(sendbuf,sendcount,datatype,succ,sendtag,
    		 recvbuf,recvcount,datatype,pred,MPI_ANY_TAG,
    		 comm,&status);
    if (pred!=MPI_PROC_NULL) {
      if (!isroot&&status.MPI_TAG==BCSTTAG) isroot = 1; // now root
    }
  }

  free(nullbuf);

  return MPI_SUCCESS;
}

int comp_base(const void *a, const void *b)
{
  const base_t *da = (const base_t*)a;
  const base_t *db = (const base_t*)b;

  return (*da > *db) - (*da < *db);
}

#define MEDIAN3(_a0,_a1,_a2) \
((_a0 > _a1) ? (_a1 > _a2) ? _a1 : (_a0 > _a2) ? _a2 : _a0 : (_a1 > _a2) ? (_a0 > _a2) ? _a0 : _a2 : _a1)

#define MEDIAN5(_a0,_a1,_a2,_a4,_a5)

inline base_t choose(base_t *a, int n)
{
  base_t pivot;

#ifdef USEINTERPOLATION
  pivot = a[random()%n];
#else
  pivot = MEDIAN3(a[random()%n],a[random()%n],a[random()%n]);
#endif
  
  return pivot;
}

inline int ranindicator(const int p)
{
  const int pool = RANGE;
  const int range = pool*pool;
  int lottery = 1+p/pool;
  int indicator;

  if (random()%lottery==0) {
    indicator = 1+random()%range;
  } else indicator = 1;
  
  return indicator;
}

inline void select_minmax(base_t *a, int n, base_t minmaxelt[2])
{
  int i;
  
  if (n>0) {
    if (a[0]<a[n-1]) {
      minmaxelt[0] = a[0];
      minmaxelt[1] = a[n-1];
    } else {
      minmaxelt[0] = a[n-1];
      minmaxelt[1] = a[0];
    }
    for (i=1; i<=MMS&&i<n; i++) {
      if (a[i]<minmaxelt[0]) minmaxelt[0] = a[i]; else {
	if (a[i]>minmaxelt[1]) minmaxelt[1] = a[i];
      }
    }
    for (i=n-2; i>=n-2-MMS&&i>0; i--) {
      if (a[i]<minmaxelt[0]) minmaxelt[0] = a[i]; else {
	if (a[i]>minmaxelt[1]) minmaxelt[1] = a[i];
      }
    }

    assert(minmaxelt[0]<=minmaxelt[1]);
  } else {
    minmaxelt[0] = MAXBASE;
    minmaxelt[1] = -MAXBASE;
  }
}
	   
// for improved load balancing, exploit partition into <, =, > better
// Bentley, McIlroy '93

#define SWAP(_a,_i,_j) \
{register base_t tmp = _a[_i]; _a[_i] = _a[_j]; _a[_j] = tmp; }

int partition(base_t a[], const int n, const base_t pivot)
{
  register int i, j;

  if (n==0) return 0;
  if (n==1) {
    if (a[0]<=pivot) return 1; else return 0;
  }
 
#ifdef USESENTINEL
  base_t a0, an;
  int swap0, swapn;

  i = 0; j = n-1;
  a0 = a[0]; an = a[j];
  if (a0>=pivot) {
    a[0] = pivot; swap0 = 1;
  } else swap0 = 0;
  if (an<=pivot) {
    a[j] = pivot; swapn = 1;
  } else swapn = 0;
#else
  i = -1; j = n;
#endif

  // Invariant: a[0...i]<=pivot, a[j...n-1]>=pivot

  do {
#ifdef USESENTINEL
    while (a[++i]<pivot);
    while (a[--j]>pivot);
#else
    while (i<j&&a[++i]<pivot); // could be improved by double sentinel?
    while (i<j&&a[--j]>pivot);
#endif
    if (i>=j) break; // invariant violated
    SWAP(a,i,j);     // restablish
  } while (1);

#ifdef USESENTINEL
  /* Single sentinel
  if (a0<=pivot) a[0] = a0; else {
    i--;
    assert(a[i]<=pivot);
    a[0] = a[i]; a[i] = a0;
  }
  */
  if (swap0) {
    if (swapn) {
      a[n-1] = a0; a[0] = an;
    } else {
      i--;
      assert(a[i]<=pivot);
      a[0] = a[i]; a[i] = a0;
    }
  } else if (swapn) {
    assert(a[i]>=pivot);
    a[n-1] = a[i]; a[i] = an; i++;
  }
#endif
  /*
  if (i<n) {
    while (i>n/2&&a[i]==pivot) i--;
  }
  if (i>=0) {
    while (i<n/2&&a[i]==pivot) i++;
  }
  */
  return i;
}

int split(const base_t a[], const int n, const base_t pivot)
{
  // by binary search
  register int i, j, k;

  if (n==0) return 0;
  if (a[0]>pivot) return 0;
  if (a[n-1]<pivot) return n;
  i = 0; j = n;
  while (i<j-1) {
    k = (i+j)/2;
    if (a[k]==pivot) {
      if (k<n/2) i = k; else j = k; // split towards middle 
    } else if (a[k]<pivot) i = k; else j = k;
  }
  assert(a[j-1]<=pivot&&(j==n||a[j]>=pivot));

  return j;
}

#ifdef STANDARDMERGE
int merge(const base_t a[], const int n, const base_t b[], const int m, 
	  base_t c[])
{
  register int i, j, k;

  i = 0; j = 0; k = 0;
  while (i<n&&j<m) {
    c[k++] = (a[i]<=b[j]) ? a[i++] : b[j++];
  }

  while (i<n) c[k++] = a[i++]; 
  while (j<m) c[k++] = b[j++]; 
    
  return k;
}
#else
int merge(const base_t a[], const int n, const base_t b[], const int m, 
	  base_t c[])
{
  register int i, j, k;

  register base_t x;

  i = 0; j = 0; k = 0;
  if (n>0&&m>0) {
    if (a[i]<=b[j]) {
      do {
	x = b[j];
	do { 
	  c[k++] = a[i++];
	} while (i<n&&a[i]<=x);
	if (i<n) x = a[i]; else x = b[m-1];
	do {
	  c[k++] = b[j++];
	} while (j<m&&b[j]<=x);
      } while (i<n&&j<m);
    } else {
      do {
	x = a[i];
	do { 
	  c[k++] = b[j++];
	} while (j<m&&b[j]<x);
	if (j<m) x = b[j]; else x = a[n-1];
	do {
	  c[k++] = a[i++];
	} while (i<n&&a[i]<=x);
      } while (i<n&&j<m);
    }
  }
  while (i<n) c[k++] = a[i++]; 
  while (j<m) c[k++] = b[j++]; 

  return k;
}
#endif

// n-way heap-based merge, heap-oredred array, heap stores indices
inline int siftup(int h, int heap[], const base_t a[], const int as[])
{
  while (h>1&&a[as[heap[h/2]]]>a[as[heap[h]]]) {
    SWAP(heap,h,h/2);
    h = h/2;
  }
}

inline int siftdown(int k, 
		    const int h, int heap[], const base_t a[], const int as[])
{
  while (2*k<=h) {
    int j = 2*k;
    if (j<h&&a[as[heap[j]]]>a[as[heap[j+1]]]) j++;
    if (a[as[heap[k]]]<=a[as[heap[j]]]) break;
    SWAP(heap,j,k);
    k = j;
  }
}

inline int heapdelmin(const base_t a[], const int as[], int *h, int heap[])
{
  int i;

  assert(*h>0);
  i = heap[1]; // assume non-empty
  heap[1] = heap[*h];
  (*h)--;
  siftdown(1,*h,heap,a,as);

  return i;
}

inline int heapins(int i, const base_t a[], const int as[], int h, int heap[])
{
  heap[++h] = i;
  siftup(h,heap,a,as);

  return h;
}

void heapify(const base_t a[], const int as[], int h, int heap[])
{
  register int i;

  for (i=h/2+1; i>0; i--) {
    siftdown(i,h,heap,a,as);
  }
}

int multimerge(const base_t a[], const int n, int as[], const int au[],
	       base_t b[])
{
  register int i, j, k;
  int *heap, h;

  heap = (int*)malloc((n+1)*sizeof(int));
  h = 0;

#ifdef STANDARDMERGE
  for (i=0; i<n; i++) {
    if (as[i]<au[i]) h = heapins(i,a,as,h,heap);
  }

  k = 0;
  while (h>0) {
    i = heapdelmin(a,as,&h,heap);
    b[k++] = a[as[i]]; as[i]++;
    if (as[i]<au[i]) h = heapins(i,a,as,h,heap);
  }
#else
  for (i=0; i<n; i++) {
    if (as[i]<au[i]) heap[++h] = i;
  }
  heapify(a,as,h,heap); // O(h)

  k = 0;
  if (h>0) {
    i = heapdelmin(a,as,&h,heap);

    do {
      if (h>0) j = heapdelmin(a,as,&h,heap); else break;
      register base_t x = a[as[j]];
      // postpone going to heap as far as possible
      while (as[i]<au[i]&&a[as[i]]<=x) b[k++] = a[as[i]++];
      if (as[i]<au[i]) h = heapins(i,a,as,h,heap);
      i = j;
    } while (1);
    while (as[i]<au[i]) b[k++] = a[as[i]++];
  }
#endif

  free(heap);

  return k;
}

// standard hypercube Quicksort
// must be called with *m==0 to prevent a from being freed
#ifdef COMSORT
base_t *par_qsort(base_t *a, int n, int *m, MPI_Comm comm, int ncut, int nlog)
#else
base_t *par_qsort(base_t *a, int n, int *m, MPI_Comm comm)
#endif
{
#ifdef USEINTERPOLATION
  base_t pivot[2], localpivot[2];
#else
  baseint_t pivot[K], localpivot[K];
#endif
  base_t *b, *c;
  int nn, ns, nl;

  MPI_Comm newcomm;

  int rank, size;
  int half;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  assert(((size-1)&size)==0); // size must be power of two

  if (size==1) {
    *m = n;
    qsort(a,n,sizeof(base_t),comp_base);

    return a;
  }

#ifdef COMSORT
  if (size>2&&ncut>PARFACT*size/nlog) return tuw_qsort(a,n,m,comm);
#endif
      
#ifdef USEINTERPOLATION
  select_minmax(a,n,localpivot);
  localpivot[0] = -localpivot[0];
  MPI_Allreduce(&localpivot,&pivot,2,BASE_T,MPI_MAX,comm);
  pivot[0] = -pivot[0];
  pivot[0] = (pivot[0]+pivot[1])/2;
  nn = partition(a,n,pivot[0]);
#else
  if (n>0) {
    int k;
    for (k=0; k<K; k++) {
      localpivot[k].val = choose(a,n);
      localpivot[k].loc = ranindicator(size);
    }
  } else {
    int k;
    for (k=0; k<K; k++) {
      localpivot[k].val = 0;
      localpivot[k].loc = 0;
    }
  }
  MPI_Allreduce(localpivot,pivot,K,BASEINT_T,TUW_CHOOSE,comm);
  if (n>0) {
#if K==3
    base_t globalpivot = MEDIAN3(pivot[0].val,pivot[1].val,pivot[2].val);
#else
    base_t globalpivot = pivot[0].val;
#endif
    nn = partition(a,n,globalpivot);
  } else {
    nn = 0;
  }
#endif
  
  half = size/2;
  if (rank<half) {
    // smaller half (keep and get smaller elements)
    nl = n-nn;
    MPI_Sendrecv(&nl,1,MPI_INT,rank+half,QTAG,
		 &ns,1,MPI_INT,rank+half,QTAG,comm,MPI_STATUS_IGNORE);
    n = nn+ns;
    b = (base_t*)malloc(n*sizeof(base_t));
    assert(n==0||b!=NULL);
#ifdef USEBLOCKING
    MPI_Sendrecv(a+nn,nl,BASE_T,rank+half,QTAG,
		 b+nn,ns,BASE_T,rank+half,QTAG,comm,MPI_STATUS_IGNORE);
    memcpy(b,a,nn*sizeof(base_t)); // use loop for element copy (or avoid)
#else
    MPI_Request request[4];
    MPI_Isend(a+nn,nl,BASE_T,rank+half,QTAG,comm,&request[0]);
    MPI_Irecv(b+nn,ns,BASE_T,rank+half,QTAG,comm,&request[1]);
#ifdef USEMEMCPY
    memcpy(b,a,nn*sizeof(base_t)); // use loop for element copy (or avoid)
    MPI_Waitall(2,request,MPI_STATUSES_IGNORE);
#else
    MPI_Isend(a,nn,BASE_T,0,QTAG,MPI_COMM_SELF,&request[2]);
    MPI_Irecv(b,nn,BASE_T,0,QTAG,MPI_COMM_SELF,&request[3]);
    MPI_Waitall(4,request,MPI_STATUSES_IGNORE);
#endif
#endif
  } else {
    // larger half (keep and get larger elements)
    ns = nn;
    MPI_Sendrecv(&ns,1,MPI_INT,rank-half,QTAG,
		 &nl,1,MPI_INT,rank-half,QTAG,comm,MPI_STATUS_IGNORE);
    n = n-nn+nl;
    b = (base_t*)malloc(n*sizeof(base_t));
    assert(n==0||b!=NULL);
#ifdef USEBLOCKING
    MPI_Sendrecv(a,ns,BASE_T,rank-half,QTAG,
		 b,nl,BASE_T,rank-half,QTAG,comm,MPI_STATUS_IGNORE);
    memcpy(b+nl,a+ns,(n-nl)*sizeof(base_t)); // use loop for element copy (or avoid)
#else
    MPI_Request request[4];
    MPI_Isend(a,ns,BASE_T,rank-half,QTAG,comm,&request[0]);
    MPI_Irecv(b,nl,BASE_T,rank-half,QTAG,comm,&request[1]);
#ifdef USEMEMCPY
    memcpy(b+nl,a+ns,(n-nl)*sizeof(base_t)); // use loop for element copy (or avoid)
    MPI_Waitall(2,request,MPI_STATUSES_IGNORE);
#else
    MPI_Isend(a+ns,(n-nl),BASE_T,0,QTAG,MPI_COMM_SELF,&request[2]);
    MPI_Irecv(b+nl,(n-nl),BASE_T,0,QTAG,MPI_COMM_SELF,&request[3]);
    MPI_Waitall(4,request,MPI_STATUSES_IGNORE);
#endif
#endif
  }

  if (*m>0) free(a); 
  *m = n;

#ifdef USEGROUPS
  MPI_Group group, newgroup;
  MPI_Comm_group(comm,&group); // get group
  if (rank<half) {
    int range[1][3] = {0,half-1,1};
    MPI_Group_range_incl(group,1,range,&newgroup);
  } else {
    int range[1][3] = {half,size-1,1};
    MPI_Group_range_incl(group,1,range,&newgroup);
  }
#ifdef USEGROUPCREATE
  MPI_Comm_create_group(comm,newgroup,QTAG,&newcomm);
#else
  MPI_Comm_create(comm,newgroup,&newcomm);
#endif
  MPI_Group_free(&group);
  MPI_Group_free(&newgroup);
#else
  int color = (rank<half) ? 0 : 1;
  MPI_Comm_split(comm,color,0,&newcomm);
#endif

#ifdef COMSORT
  c = par_qsort(b,n,m,newcomm,ncut,nlog-1);
#else
  c = par_qsort(b,n,m,newcomm);
#endif

  MPI_Comm_free(&newcomm);

  return c;
}

// HyperQuicksort, sort-first, then recursive merge exchange
#ifdef COMSORT
base_t *mergepartition(base_t *a, int n, int *m, MPI_Comm comm,
		       int ncut, int nlog)
#else
base_t *mergepartition(base_t *a, int n, int *m, MPI_Comm comm)
#endif
{
#ifdef USEINTERPOLATION
  base_t pivot[2], localpivot[2];
#else
  baseint_t pivot[K], localpivot[K];
#endif
  base_t *b, *c;
  int nn, ns, nl;

  MPI_Comm newcomm;

  int rank, size;
  int half;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  assert(((size-1)&size)==0); // size must be power of two

  if (size==1) {
    *m = n; return a;
  }

#ifdef COMSORT
  if (size>2&&ncut>HYPFACT*size/nlog) return wut_merge(a,n,m,comm);
#endif

#ifdef USEINTERPOLATION
  select_minmax(a,n,localpivot);
  localpivot[0] = -localpivot[0];
  MPI_Allreduce(&localpivot,&pivot,2,BASE_T,MPI_MAX,comm);
  pivot[0] = -pivot[0];
  pivot[0] = (pivot[0]+pivot[1])/2;
  nn = split(a,n,pivot[0]); // partition of sorted array is easy
#else
  if (n>0) {
    int k;
    for (k=0; k<K; k++) {
      localpivot[k].val = a[n/2]; 
      localpivot[k].loc = ranindicator(size);
    }
  } else {
    int k;
    for (k=0; k<K; k++) {
      localpivot[k].val = 0;
      localpivot[k].loc = 0;
    }
  }
  MPI_Allreduce(localpivot,pivot,K,BASEINT_T,TUW_CHOOSE,comm);
  if (n>0) {
#if K==3
    base_t globalpivot = MEDIAN3(pivot[0].val,pivot[1].val,pivot[2].val);
#else
    base_t globalpivot = pivot[0].val;
#endif
    nn = split(a,n,globalpivot); // partition of sorted array is easy
  } else {
    nn = 0;
  }
#endif
  assert(0<=nn&&nn<=n);
  
  half = size/2;
  if (rank<half) {
    // smaller half (keep and get smaller elements)
    nl = n-nn;
    MPI_Sendrecv(&nl,1,MPI_INT,rank+half,QTAG,
		 &ns,1,MPI_INT,rank+half,QTAG,comm,MPI_STATUS_IGNORE);
    n = nn+ns;
    b = (base_t*)malloc(ns*sizeof(base_t));
    c = (base_t*)malloc(n*sizeof(base_t));
    assert(n==0||b!=NULL);
    assert(n==0||c!=NULL);

    MPI_Sendrecv(a+nn,nl,BASE_T,rank+half,QTAG,
		 b,ns,BASE_T,rank+half,QTAG,comm,MPI_STATUS_IGNORE);
    int nnn = merge(a,nn,b,ns,c);
    assert(nnn==n);
  } else {
    // larger half (keep and get larger elements)
    ns = nn;
    MPI_Sendrecv(&ns,1,MPI_INT,rank-half,QTAG,
		 &nl,1,MPI_INT,rank-half,QTAG,comm,MPI_STATUS_IGNORE);
    n = n-nn+nl;
    b = (base_t*)malloc(nl*sizeof(base_t));
    c = (base_t*)malloc(n*sizeof(base_t));
    assert(n==0||b!=NULL);
    assert(n==0||c!=NULL);

    MPI_Sendrecv(a,ns,BASE_T,rank-half,QTAG,
		 b,nl,BASE_T,rank-half,QTAG,comm,MPI_STATUS_IGNORE);
    int nnn = merge(a+nn,n-nl,b,nl,c);
    assert(nnn==n);
  }
  free(b);

  if (*m>0) free(a);
  *m = n;

#ifdef USEGROUPS
  MPI_Group group, newgroup;
  MPI_Comm_group(comm,&group); // get group
  if (rank<half) {
    int range[1][3] = {0,half-1,1};
    MPI_Group_range_incl(group,1,range,&newgroup);
  } else {
    int range[1][3] = {half,size-1,1};
    MPI_Group_range_incl(group,1,range,&newgroup);
  }
#ifdef USEGROUPCREATE
  MPI_Comm_create_group(comm,newgroup,QTAG,&newcomm);
#else
  MPI_Comm_create(comm,newgroup,&newcomm);
#endif
  MPI_Group_free(&group);
  MPI_Group_free(&newgroup);
#else
  int color = (rank<half) ? 0 : 1;
  MPI_Comm_split(comm,color,0,&newcomm);
#endif

#ifdef COMSORT
  base_t *d = mergepartition(c,n,m,newcomm,ncut,nlog-1);
#else
  base_t *d = mergepartition(c,n,m,newcomm);
#endif
  
  MPI_Comm_free(&newcomm);

  return d;
}

#ifdef COMSORT
base_t *hyp_qsort(base_t *a, const int n, int *m, MPI_Comm comm,
		  int ncut, int nlog)
#else
base_t *hyp_qsort(base_t *a, const int n, int *m, MPI_Comm comm)
#endif
{
  qsort(a,n,sizeof(base_t),comp_base);
#ifdef COMSORT
  return mergepartition(a,n,m,comm,ncut,nlog);
#else
  return mergepartition(a,n,m,comm);
#endif
}

// partition without element exchange
base_t *tuw_qsort(base_t *a, const int n, int *m, MPI_Comm comm)
{
  int i, j, k;

  int rank, size;

#ifdef USEINTERPOLATION
  base_t *pivot, *localpivot;
#else
  baseint_t *pivot, *localpivot;
#endif
  int *adisp; // start displacements
  int *an;
  int *bdisp; // start displacements
  int *bn;

  int worka = (*m>0); // a is an intermediate work array, free
  base_t *b;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  assert(((size-1)&size)==0); // size must be power of two

#ifdef USEINTERPOLATION
  pivot = (base_t*)malloc(4*size*sizeof(base_t));
  assert(pivot!=NULL);
  localpivot = pivot+2*size;
#else
  pivot = (baseint_t*)malloc(2*K*size*sizeof(baseint_t));
  assert(pivot!=NULL);
  localpivot = pivot+K*size;
#endif
  adisp = (int*)malloc(4*size*sizeof(int));
  assert(adisp!=NULL);
  an = adisp+size;
  bdisp = an+size;
  bn = bdisp+size;

  an[0] = n;
  adisp[0] = 0;
  
  for (k=size; k>1; k>>=1) {
    // select local pivots
#ifdef USEINTERPOLATION
    for (i=0, j=0; i<size; i+=k,j+=2) {
      select_minmax(a+adisp[i],an[i],localpivot+j);
      localpivot[j] = -localpivot[j];
    }
#else
    for (i=0, j=0; i<size; i+=k,j+=K) {
      if (an[i]>0) {
	int kk;
	for (kk=0; kk<K; kk++) {
	  localpivot[j+kk].val = choose(a+adisp[i],an[i]);
	  localpivot[j+kk].loc = ranindicator(size);
	}
      } else {
	int kk;
	for (kk=0; kk<K; kk++) {
	  localpivot[j+kk].val = 0;
	  localpivot[j+kk].loc = 0;
	}
      }
    }
#endif

    // compute global pivots
#ifdef USEINTERPOLATION
    MPI_Allreduce(localpivot,pivot,j,BASE_T,MPI_MAX,comm);
#else
    MPI_Allreduce(localpivot,pivot,j,BASEINT_T,TUW_CHOOSE,comm);
#endif

    // distributed partition
    int checkn = 0;
#ifdef USEINTERPOLATION
    for (i=0, j=0; i<size; i+=k,j+=2) {
      int inext = i+(k>>1);
      if (an[i]>0) {
	pivot[j] = -pivot[j];
	pivot[j] = (pivot[j]+pivot[j+1])/2; 
	adisp[inext]  = partition(a+adisp[i],an[i],pivot[j]);
	an[inext]     = an[i]-adisp[inext];
	an[i]         = adisp[inext];
	adisp[inext] += adisp[i];
      } else {
	an[inext] = 0;
	adisp[inext] = adisp[i];
      }
      checkn += an[i]+an[inext];
    }
#else
    for (i=0, j=0; i<size; i+=k,j+=K) {
      int inext = i+(k>>1);
      if (an[i]>0) {
	base_t globalpivot;
#if K==3
	globalpivot = MEDIAN3(pivot[j].val,pivot[j+1].val,pivot[j+2].val);
#else
	globalpivot = pivot[j].val;
#endif
	adisp[inext]  = partition(a+adisp[i],an[i],globalpivot);
	an[inext]     = an[i]-adisp[inext];
	an[i]         = adisp[inext];
	adisp[inext] += adisp[i];
      } else {
	an[inext] = 0;
	adisp[inext] = adisp[i];
      }
      checkn += an[i]+an[inext];
    }
#endif
    assert(checkn==n);
  }

  MPI_Alltoall(an,1,MPI_INT,bn,1,MPI_INT,comm);
  bdisp[0] = 0;
  for (i=1; i<size; i++) {
    bdisp[i] = bdisp[i-1]+bn[i-1];
  }
  *m = bdisp[size-1]+bn[size-1];
  b = (base_t*)malloc(*m*sizeof(base_t));
  assert(b!=NULL);

  MPI_Alltoallv(a,an,adisp,BASE_T,b,bn,bdisp,BASE_T,comm);
  if (worka) free(a);

  // Quicksort locally
  qsort(b,*m,sizeof(base_t),comp_base);
  
  free(pivot);
  free(adisp);
  
  return b;
}

// mergepartition without element exchange
base_t *wut_merge(base_t *a, const int n, int *m, MPI_Comm comm)
{
  int i, j, k;

  int rank, size;

#ifdef USEINTERPOLATION
  base_t *pivot, *localpivot;
#else
  baseint_t *pivot, *localpivot;
#endif
  int *adisp; // start displacements
  int *an;
  int *bdisp; // start displacements
  int *bn;

  int worka = (*m>0); // a is an intermediate work array, free
  base_t *b, *c;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  assert(((size-1)&size)==0); // size must be power of two

#ifdef USEINTERPOLATION
  pivot = (base_t*)malloc(4*size*sizeof(base_t));
  assert(pivot!=NULL);
  localpivot = pivot+2*size;
#else
  pivot = (baseint_t*)malloc(2*K*size*sizeof(baseint_t));
  assert(pivot!=NULL);
  localpivot = pivot+K*size;
#endif
  adisp = (int*)malloc(4*size*sizeof(int));
  assert(adisp!=NULL);
  an = adisp+size;
  bdisp = an+size;
  bn = bdisp+size;

  an[0] = n;
  adisp[0] = 0;

  for (k=size; k>1; k>>=1) {
    // select local pivots
#ifdef USEINTERPOLATION
    for (i=0, j=0; i<size; i+=k,j+=2) {
      select_minmax(a+adisp[i],an[i],localpivot+j);
      localpivot[j] = -localpivot[j];
    }
#else
    for (i=0, j=0; i<size; i+=k,j+=K) {
      if (an[i]>0) {
	int kk;
	for (kk=0; kk<K; kk++) {
	  localpivot[j+kk].val = a[adisp[i]+an[i]/2];
	  localpivot[j+kk].loc = ranindicator(size);
	}
      } else {
	int kk;
	for (kk=0; kk<K; kk++) {
	  localpivot[j+kk].val = 0;
	  localpivot[j+kk].loc = 0;
	}
      }
    }
#endif

    // compute global pivots (idea: use average)
#ifdef USEINTERPOLATION
    MPI_Allreduce(localpivot,pivot,j,BASE_T,MPI_MAX,comm);
#else
    MPI_Allreduce(localpivot,pivot,j,BASEINT_T,TUW_CHOOSE,comm);
#endif

    // distributed partition
    int checkn = 0;
#ifdef USEINTERPOLATION
    for (i=0, j=0; i<size; i+=k,j+=2) {
      int inext = i+(k>>1);
      if (an[i]>0) {
	pivot[j] = -pivot[j];
	pivot[j] = (pivot[j]+pivot[j+1])/2;
	adisp[inext]  = split(a+adisp[i],an[i],pivot[j]);
	an[inext]     = an[i]-adisp[inext];
	an[i]         = adisp[inext];
	adisp[inext] += adisp[i];
      } else {
	an[inext] = 0;
	adisp[inext] = adisp[i];
      }
      checkn += an[i]+an[inext];
    }
#else
    for (i=0, j=0; i<size; i+=k,j+=K) {
      int inext = i+(k>>1);
      if (an[i]>0) {
	base_t globalpivot;
#if K==3
	globalpivot = MEDIAN3(pivot[j].val,pivot[j+1].val,pivot[j+2].val);
#else
	globalpivot = pivot[j].val;
#endif
	adisp[inext]  = split(a+adisp[i],an[i],globalpivot);
	an[inext]     = an[i]-adisp[inext];
	an[i]         = adisp[inext];
	adisp[inext] += adisp[i];
      } else {
	an[inext] = 0;
	adisp[inext] = adisp[i];
      }
      checkn += an[i]+an[inext];
    }
#endif
    assert(checkn==n);
  }

  MPI_Alltoall(an,1,MPI_INT,bn,1,MPI_INT,comm);
  bdisp[0] = 0;
  for (i=1; i<size; i++) {
    bdisp[i] = bdisp[i-1]+bn[i-1];
  }
  *m = bdisp[size-1]+bn[size-1];
  b = (base_t*)malloc(*m*sizeof(base_t));
  assert(b!=NULL);

  MPI_Alltoallv(a,an,adisp,BASE_T,b,bn,bdisp,BASE_T,comm);
  if (worka) free(a);
  
  // p-way merge
  c = (base_t*)malloc(*m*sizeof(base_t));
  assert(c!=NULL);

  for (i=0; i<size; i++) {
    int nn = bn[i];
    bn[i] = bdisp[i]; bdisp[i] += nn;
  }
  
  int checkn = multimerge(b,size,bn,bdisp,c);
  assert(checkn==*m);

  free(b);

  free(pivot);
  free(adisp);
  
  return c;
}

base_t *wut_qsort(base_t *a, const int n, int *m, MPI_Comm comm)
{
  // Quicksort locally
  qsort(a,n,sizeof(base_t),comp_base);

  return wut_merge(a,n,m,comm);
}

// Iterative Quicksort (no communicator split)
base_t *ite_qsort(base_t *a, int n, int *m, MPI_Comm comm)
{
#ifdef USEINTERPOLATION
  basebase_t pivot, segmentpivot;
#else
  baseint_t pivot[K], segmentpivot[K];
#endif
  base_t globalpivot;

  base_t *b;
  int nn, ns, nl;

  MPI_Comm newcomm;

  int rank, size;

  int k, kk;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  assert(((size-1)&size)==0); // size must be power of two
  
  *m = 0;
  for (k=size; k>1; k = kk) {
    if (n>0) {
#ifdef USEINTERPOLATION
      pivot.val = choose(a,n);
      pivot.loc = 1;
#else
      int j;
      for (j=0; j<K; j++) {
	pivot[j].val = choose(a,n);
	pivot[j].loc = ranindicator(size)+(rank/k)*(RANGE*RANGE);
      }
#endif
    } else {
#ifdef USEINTERPOLATION
      pivot.val = 0;
      pivot.loc = 0;
#else
      int j;
      for (j=0; j<K; j++) {
	pivot[j].val = 0;
	pivot[j].loc = 0;
      }
#endif
    }
#ifdef USEINTERPOLATION
    // Segmented prefix-sums like computation of pivot average
    MPI_Scan(&pivot,&segmentpivot,2,BASE_T,MPI_SUM,comm);
#else
    MPI_Scan(pivot,segmentpivot,K,BASEINT_T,TUW_CHOOSE,comm);
#endif
    if ((rank&(k-1))==k-1) {
#ifdef USEINTERPOLATION
      MPI_Send(&segmentpivot,2,BASE_T,rank-(k-1),QTAG,comm);
#else
      MPI_Send(&segmentpivot,K,BASEINT_T,rank-(k-1),QTAG,comm);
#endif
    } else if ((rank&(k-1))==0x0) {
#ifdef USEINTERPOLATION
      if (rank!=0) {
	segmentpivot.val -= pivot.val;
	segmentpivot.loc -= pivot.loc;
      } else {
	segmentpivot.val = 0;
	segmentpivot.loc = 0;
      }
      MPI_Recv(&pivot,2,BASE_T,rank+(k-1),QTAG,comm,MPI_STATUS_IGNORE);
      pivot.val -= segmentpivot.val;
      pivot.loc -= segmentpivot.loc;
      if (pivot.loc>=1) pivot.val /= pivot.loc;
      globalpivot = pivot.val;
#else
      MPI_Recv(&pivot,K,BASEINT_T,rank+(k-1),QTAG,comm,MPI_STATUS_IGNORE);
#if K==3
      globalpivot = MEDIAN3(pivot[0].val,pivot[1].val,pivot[2].val);
#else
      globalpivot = pivot[0].val;
#endif
      
#endif
    }
#ifdef USEBUILTIN
    baseint_t segpivot;
    segpivot.val = globalpivot;
    segpivot.loc = ((rank&(k-1))==0x0) ? 1 : 0; // mark segment starts
    MPI_Scan(MPI_IN_PLACE,&segpivot,1,BASEINT_T,TUW_SEGREP,comm);
    globalpivot = segpivot.val;
#else
    int root = ((rank&(k-1))==0x0); // segment starts
    bcast_segment(&globalpivot,1,BASE_T,root,comm);
#endif

    if (n>0) {
      nn = partition(a,n,globalpivot);
    } else {
      nn = 0;
    }

    kk = k>>1;
    if ((rank&kk)!=kk) {
      // smaller half (keep and get smaller elements)
      nl = n-nn;
      MPI_Sendrecv(&nl,1,MPI_INT,rank+kk,QTAG,
		   &ns,1,MPI_INT,rank+kk,QTAG,comm,MPI_STATUS_IGNORE);
      
      n = nn+ns;
      b = (base_t*)malloc(n*sizeof(base_t));
#ifdef USEBLOCKING
      MPI_Sendrecv(a+nn,nl,BASE_T,rank+kk,QTAG,
		   b+nn,ns,BASE_T,rank+kk,QTAG,comm,MPI_STATUS_IGNORE);
      memcpy(b,a,nn*sizeof(base_t)); // use loop for element copy (or avoid)
#else
      MPI_Request request[4];
      MPI_Isend(a+nn,nl,BASE_T,rank+kk,QTAG,comm,&request[0]);
      MPI_Irecv(b+nn,ns,BASE_T,rank+kk,QTAG,comm,&request[1]);
#ifdef USEMEMCPY
      memcpy(b,a,nn*sizeof(base_t)); // use loop for element copy (or avoid)
      MPI_Waitall(2,request,MPI_STATUSES_IGNORE);
#else
      MPI_Isend(a,nn,BASE_T,0,QTAG,MPI_COMM_SELF,&request[2]); 
      MPI_Irecv(b,nn,BASE_T,0,QTAG,MPI_COMM_SELF,&request[3]);
      MPI_Waitall(4,request,MPI_STATUSES_IGNORE);
#endif
#endif
    } else {
      // larger k (keep and get larger elements)
      ns = nn;
      MPI_Sendrecv(&ns,1,MPI_INT,rank-kk,QTAG,
		   &nl,1,MPI_INT,rank-kk,QTAG,comm,MPI_STATUS_IGNORE);
      n = n-nn+nl;
      b = (base_t*)malloc(n*sizeof(base_t));
#ifdef USEBLOCKING
      MPI_Sendrecv(a,ns,BASE_T,rank-kk,QTAG,
		   b,nl,BASE_T,rank-kk,QTAG,comm,MPI_STATUS_IGNORE);
      memcpy(b+nl,a+ns,(n-nl)*sizeof(base_t)); // use loop for element copy (or avoid)
#else
      MPI_Request request[4];
      MPI_Isend(a,ns,BASE_T,rank-kk,QTAG,comm,&request[0]);
      MPI_Irecv(b,nl,BASE_T,rank-kk,QTAG,comm,&request[1]);
#ifdef USEMEMCPY
      memcpy(b+nl,a+ns,(n-nl)*sizeof(base_t)); // use loop for element copy (or avoid)
      MPI_Waitall(2,request,MPI_STATUSES_IGNORE);
#else
      MPI_Isend(a+ns,(n-nl),BASE_T,0,QTAG,MPI_COMM_SELF,&request[2]);
      MPI_Irecv(b+nl,(n-nl),BASE_T,0,QTAG,MPI_COMM_SELF,&request[3]);
      MPI_Waitall(4,request,MPI_STATUSES_IGNORE);
#endif
#endif
    }
    if (*m>0) free(a);
    a = b;
    *m = n;
  }
  *m = n;

  qsort(a,n,sizeof(base_t),comp_base);

  return a;
}

// predicate
#define SORTTAG 1111
int sorted(base_t a[], int n, MPI_Comm comm)
{
  int rank, size;
  int good;
  
  int i;

  base_t ap;
  
  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);
  
  for (i=0; i<n-1; i++) if (a[i]>a[i+1]) break;
  good = (i==n-1);

  if (size>1) {
    if (rank==0) {
      MPI_Send(&a[n-1],1,BASE_T,rank+1,SORTTAG,comm);
    } else if (rank==size-1) {
      MPI_Recv(&ap,1,BASE_T,rank-1,SORTTAG,comm,MPI_STATUS_IGNORE);
      good = good&&(ap<=a[0]);
    } else {
      MPI_Sendrecv(&a[n-1],1,BASE_T,rank+1,SORTTAG,
		   &ap,1,BASE_T,rank-1,SORTTAG,comm,MPI_STATUS_IGNORE);
      good = good&&(ap<=a[0]);
    }
  }
  
  return good;
}

int fisheryates(base_t a[], int n)
{
  int i, j;
  base_t aa;
  
  for (i=0; i<n-1; i++) {
    j = i+random()%(n-i);
    aa = a[i]; a[i] = a[j]; a[j] = aa;
  }
}

/* Simplified Sanders 1998 */
int ranpermutation(base_t a[], int n, MPI_Comm comm)
{
  int rank, size;
  int i;
  base_t *b;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  assert(n%size==0);

  b = (base_t*)malloc(n*sizeof(base_t));

  for (i=0; i<n; i++) b[i] = a[i];
  //fisheryates(a,n); // permute locally

  MPI_Alltoall(b,n/size,BASE_T,
	       a,n/size,BASE_T,comm);

  free(b);

  fisheryates(a,n); // permute locally
}

base_t *generate(int n, int inputtype, long range, long seed)
{
  int i;
  base_t *a;

  a = (base_t*)malloc(n*sizeof(base_t));
  assert(a!=NULL);

  srandom(seed);
  if (inputtype==ASC) {
    // ascending
    int k = 0;
    int s;
    s = random()%range;
    for (i=0; i<n; i++) {
      if (i>=s) {
	s += random()%range;
	k++;
      }
      a[i] = (base_t)k;
    }
  } else if (inputtype==PRM) {
    // permutation
    for (i=0; i<n; i++) a[i] = (base_t)i;
    fisheryates(a,n);
  } else {
    // default
    for (i=0; i<n; i++) {
      a[i] = (base_t)(random()%range);
    }
  }

  return a;
}

int main(int argc, char *argv[])
{
  int rank, size;

  int n, m;
  int np, mp;
  int i, k;
  int r, rr;

  int ncut, nlog; // Quicksort cut off for switching to exchange-free
  
  base_t *a, *b, *c, *d;
  base_t *aa;

  int *count;
  int *displ;

  long seed, range;
  int inputtype, weak;

  MPI_Init(&argc,&argv);
  MPI_Pcontrol(0);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  weak = 0; // default strong scaling

  MPI_Aint lb, extent;
  MPI_Type_get_extent(BASEINT_T,&lb,&extent);
  assert(extent==sizeof(baseint_t));

#ifdef USEINTERPOLATION
  MPI_Op_create(tuw_selsum,0,&TUW_CHOOSE);
#else
  MPI_Op_create(tuw_indmax,0,&TUW_CHOOSE);
#endif
  MPI_Op_create(tuw_segrep,0,&TUW_SEGREP);
  
  /* generate */
  if (rank==0) {
    n = 1000000;
    seed = 1;
    inputtype = UNI;
    range = n;

    for (i=1; i<argc&&argv[i][0]=='-'; i++) {
      if (argv[i][1]=='n') i++,sscanf(argv[i],"%d",&n),range = n;
      else if (argv[i][1]=='s') i++,sscanf(argv[i],"%ld",&seed);
      else if (argv[i][1]=='r') i++,sscanf(argv[i],"%ld",&range);
      else if (argv[i][1]=='A') inputtype = ASC;
      else if (argv[i][1]=='P') inputtype = PRM;
      else if (argv[i][1]=='W') weak = 1;
      else fprintf(stderr,"Unknown option %c, ignored\n",argv[i][1]);
    }
    assert(range>0);
    if (inputtype==PRM) {
      if (weak) n = size*((n+size-1)/size);
      else n = (size*size)*((n+size*size-1)/(size*size));
    }
    
    /* distribute size */
    MPI_Bcast(&n,   1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&weak,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&seed,1,MPI_LONG,0,MPI_COMM_WORLD);
    if (weak) {
      np = n;
      ncut = np;
    } else {
      np = n/size;
      ncut = np+n%size;
    }

    seed += rank;
    a = generate(n,inputtype,range,seed);

    aa = (base_t*)malloc(np*sizeof(base_t));
    assert(aa!=NULL);

    count = (int*)malloc(size*sizeof(int));
    displ = (int*)malloc(size*sizeof(int));
      
    if (!weak) {
      count[0] = np; displ[0] = 0;
      for (i=1; i<size; i++) {
	count[i] = np;
	displ[i] = displ[i-1]+count[i-1];
      }
      count[size-1] += n%size;
      MPI_Scatterv(a,count,displ,BASE_T,aa,np,BASE_T,0,MPI_COMM_WORLD);
    } else {
      MPI_Bcast(&inputtype,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&range,1,MPI_LONG,0,MPI_COMM_WORLD);

      memcpy(aa,a,np*sizeof(base_t));

      if (inputtype==PRM) {
	ranpermutation(aa,np,MPI_COMM_WORLD);
      }
    }
  } else {
    /* distribute */
    MPI_Bcast(&n,   1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&weak,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&seed,1,MPI_LONG,0,MPI_COMM_WORLD);
    if (weak) {
      np = n; 
      ncut = np;
    } else {
      np = n/size; 
      ncut = np+n%size;
      if (rank==size-1) np += n%size;
    }

    assert(np>0); // only non-trivial problems

    seed += rank;

    if (!weak) {
      aa = (base_t*)malloc(np*sizeof(base_t));
      assert(aa!=NULL);

      count = NULL;
      displ = NULL;
      
      MPI_Scatterv(a,count,displ,BASE_T,aa,np,BASE_T,0,MPI_COMM_WORLD);
    } else {
      MPI_Bcast(&inputtype,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(&range,1,MPI_LONG,0,MPI_COMM_WORLD);

      aa = generate(np,inputtype,range,seed);
      if (inputtype==PRM) {
	for (i=0; i<np; i++) aa[i] += rank*np;
	ranpermutation(aa,np,MPI_COMM_WORLD);
      }
    }
  }
  b = (base_t*)malloc(np*sizeof(base_t));
  assert(b!=NULL);

  nlog = ilog2(size)+ilog2(ncut); // log(n) for cut-off
  
  /* "The benchmark": Each implementation, one after the other */

  double start, stop;
  double seqtime, partime[REPEAT];

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Pcontrol(1);
  if (rank==0) {
    fprintf(stdout,
	    "#%s(p=%d, n=%d, np=%d, r=%ld, s=%ld, type=%d, %d repetitions):\n",
	    argv[0],size,n,np,range,seed,inputtype,REPEAT);

    c = (base_t*)malloc(n*sizeof(base_t));
    memcpy(c,a,n*sizeof(base_t));

    start = MPI_Wtime();
    qsort(a,n,sizeof(base_t),comp_base);
    stop = MPI_Wtime();
    seqtime = stop-start;
    
    for (i=1; i<SEQREPEAT; i++) {
      memcpy(a,c,n*sizeof(base_t));
      
      start = MPI_Wtime();
      qsort(a,n,sizeof(base_t),comp_base);
      stop = MPI_Wtime();
      if (stop-start<seqtime) seqtime = stop-start;
    }
    
    free(c);
  }

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef PARQSORT
  for (r=0; r<REPEAT; r++) partime[r] = 0.0;
  c = NULL;
  MPI_Barrier(MPI_COMM_WORLD);
  for (r=0,rr=0; r<REPEAT+WARMUP; r++) {
    memcpy(b,aa,np*sizeof(base_t));
    if (c!=b) free(c); // keep last
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start = MPI_Wtime();
    mp = 0;
    srandom(seed);

#ifdef COMSORT
    c = par_qsort(b,np,&mp,MPI_COMM_WORLD,ncut,nlog);
#else
    c = par_qsort(b,np,&mp,MPI_COMM_WORLD);
#endif
    
    stop = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (r<WARMUP) continue;
    partime[rr++] = stop-start;
  }
  MPI_Allreduce(MPI_IN_PLACE,partime,REPEAT,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (rank==0) {
    double paravg, parmin;

    paravg = 0.0; 
    parmin = partime[0];
    for (r=0; r<REPEAT; r++) {
      paravg += partime[r]; 
      if (partime[r]<parmin) parmin = partime[r];
    }
    paravg /= REPEAT;

    fprintf(stdout,
	    "P %d NP %d Seq %.2f PAR %.2f (avg) %.2f (min)\n",
	    size,np,seqtime*MICRO,paravg*MICRO,parmin*MICRO);
    fprintf(stdout,"#SU-PAR %.2f (avg) %.2f (min)\n",
	    seqtime/paravg,seqtime/parmin);
  }

  /* sorted? */
  sorted(c,mp,MPI_COMM_WORLD);
  
  /* compare against globally sorted */

  MPI_Gather(&mp,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD);
  if (rank==0) {
    displ[0] = 0;
    for (i=1; i<size; i++) displ[i] = displ[i-1]+count[i-1];
    m = displ[size-1]+count[size-1];

    int minload, maxload;
    minload = n; maxload = 0;
    for (i=0; i<size; i++) {
      if (count[i]<minload) minload = count[i];
      if (count[i]>maxload) maxload = count[i];
    }
    fprintf(stdout,
	    "#Load-PAR (min,exp,max): %d %d %d Imbalance %.2f\n",minload,ncut,maxload,100.0*((double)(maxload-minload)/ncut));

    if (!weak) {
      assert(m==n);
      d = (base_t*)malloc(m*sizeof(base_t));
      assert(d!=NULL);

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
      
      assert(memcmp(a,d,m*sizeof(base_t))==0);
      free(d);
    }
  } else {
    if (!weak) {
      assert(displ==NULL); assert(count==NULL); // just make sure
      d = NULL;

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
    }
  }
  if (c!=b) free(c);
#endif

#ifdef ITEQSORT
  for (r=0; r<REPEAT; r++) partime[r] = 0.0;
  c = NULL;
  MPI_Barrier(MPI_COMM_WORLD);
  for (r=0,rr=0; r<REPEAT+WARMUP; r++) {
    memcpy(b,aa,np*sizeof(base_t));
    if (c!=b) free(c); // keep last
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start = MPI_Wtime();
    mp = 0;
    srandom(seed);

    c = ite_qsort(b,np,&mp,MPI_COMM_WORLD);
    
    stop = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (r<WARMUP) continue;
    partime[rr++] = stop-start;
  }
  MPI_Allreduce(MPI_IN_PLACE,partime,REPEAT,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (rank==0) {
    double paravg, parmin;

    paravg = 0.0; 
    parmin = partime[0];
    for (r=0; r<REPEAT; r++) {
      paravg += partime[r]; 
      if (partime[r]<parmin) parmin = partime[r];
    }
    paravg /= REPEAT;

    fprintf(stdout,
	    "P %d NP %d Seq %.2f ITE %.2f (avg) %.2f (min)\n",
	    size,np,seqtime*MICRO,paravg*MICRO,parmin*MICRO);
    fprintf(stdout,"#SU-ITE %.2f (avg) %.2f (min)\n",
	    seqtime/paravg,seqtime/parmin);
  }

  /* sorted? */
  sorted(c,mp,MPI_COMM_WORLD);

  /* compare against globally sorted */

  MPI_Gather(&mp,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD);
  if (rank==0) {
    displ[0] = 0;
    for (i=1; i<size; i++) displ[i] = displ[i-1]+count[i-1];
    m = displ[size-1]+count[size-1];

    int minload, maxload;
    minload = n; maxload = 0;
    for (i=0; i<size; i++) {
      if (count[i]<minload) minload = count[i];
      if (count[i]>maxload) maxload = count[i];
    }
    fprintf(stdout,
	    "#Load-ITE (min,exp,max): %d %d %d Imbalance %.2f\n",minload,ncut,maxload,100.0*((double)(maxload-minload)/ncut));

    if (!weak) {
      assert(m==n);
      d = (base_t*)malloc(m*sizeof(base_t));
      assert(d!=NULL);

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
      
      assert(memcmp(a,d,m*sizeof(base_t))==0);
      free(d);
    }
  } else {
    if (!weak) {
      assert(displ==NULL); assert(count==NULL); // just make sure
      d = NULL;

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
    }
  }
  if (c!=b) free(c); 
#endif

#ifdef HYPQSORT
  for (r=0; r<REPEAT; r++) partime[r] = 0.0;
  c = NULL;
  MPI_Barrier(MPI_COMM_WORLD);
  for (r=0,rr=0; r<REPEAT+WARMUP; r++) {
    memcpy(b,aa,np*sizeof(base_t));
    if (c!=b) free(c); // keep last
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start = MPI_Wtime();
    mp = 0;
    srandom(seed);

#ifdef COMSORT
    c = hyp_qsort(b,np,&mp,MPI_COMM_WORLD,ncut,nlog);
#else
    c = hyp_qsort(b,np,&mp,MPI_COMM_WORLD);
#endif
    
    stop = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (r<WARMUP) continue;
    partime[rr++] = stop-start;
  }
  MPI_Allreduce(MPI_IN_PLACE,partime,REPEAT,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (rank==0) {
    double paravg, parmin;

    paravg = 0.0; 
    parmin = partime[0];
    for (r=0; r<REPEAT; r++) {
      paravg += partime[r]; 
      if (partime[r]<parmin) parmin = partime[r];
    }
    paravg /= REPEAT;

    fprintf(stdout,
	    "P %d NP %d Seq %.2f HYP %.2f (avg) %.2f (min)\n",
	    size,np,seqtime*MICRO,paravg*MICRO,parmin*MICRO);
    fprintf(stdout,"#SU-HYP %.2f (avg) %.2f (min)\n",
	    seqtime/paravg,seqtime/parmin);
  }

  /* sorted? */
  sorted(c,mp,MPI_COMM_WORLD);
  
  /* compare against globally sorted */

  MPI_Gather(&mp,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD);
  if (rank==0) {
    displ[0] = 0;
    for (i=1; i<size; i++) displ[i] = displ[i-1]+count[i-1];
    m = displ[size-1]+count[size-1];

    int minload, maxload;
    minload = n; maxload = 0;
    for (i=0; i<size; i++) {
      if (count[i]<minload) minload = count[i];
      if (count[i]>maxload) maxload = count[i];
    }
    fprintf(stdout,
	    "#Load-HYP (min,exp,max): %d %d %d Imbalance %.2f\n",minload,ncut,maxload,100.0*((double)(maxload-minload)/ncut));

    if (!weak) {
      assert(m==n);
      d = (base_t*)malloc(m*sizeof(base_t));
      assert(d!=NULL);

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
      
      assert(memcmp(a,d,m*sizeof(base_t))==0);
      free(d);
    }
  } else {
    if (!weak) {
      assert(displ==NULL); assert(count==NULL); // just make sure
      d = NULL;

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
    }
  }
  if (c!=b) free(c); 
#endif

#ifdef TUWQSORT
  for (r=0; r<REPEAT; r++) partime[r] = 0.0;
  c = NULL;
  MPI_Barrier(MPI_COMM_WORLD);
  for (r=0,rr=0; r<REPEAT+WARMUP; r++) {
    memcpy(b,aa,np*sizeof(base_t));
    if (c!=b) free(c); // keep last
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start = MPI_Wtime();
    mp = 0;
    srandom(seed);

    c = tuw_qsort(b,np,&mp,MPI_COMM_WORLD);

    stop = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (r<WARMUP) continue;
    partime[rr++] = stop-start;
  }
  MPI_Allreduce(MPI_IN_PLACE,partime,REPEAT,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (rank==0) {
    double paravg, parmin;

    paravg = 0.0; 
    parmin = partime[0];
    for (r=0; r<REPEAT; r++) {
      paravg += partime[r]; 
      if (partime[r]<parmin) parmin = partime[r];
    }
    paravg /= REPEAT;

    fprintf(stdout,
	    "P %d NP %d Seq %.2f TUW %.2f (avg) %.2f (min)\n",
	    size,np,seqtime*MICRO,paravg*MICRO,parmin*MICRO);
    fprintf(stdout,"#SU-TUW %.2f (avg) %.2f (min)\n",
	    seqtime/paravg,seqtime/parmin);
  }

  /* sorted? */
  sorted(c,mp,MPI_COMM_WORLD);
  
  /* compare against globally sorted */

  MPI_Gather(&mp,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD);
  if (rank==0) {
    displ[0] = 0;
    for (i=1; i<size; i++) displ[i] = displ[i-1]+count[i-1];
    m = displ[size-1]+count[size-1];

    int minload, maxload;
    minload = n; maxload = 0;
    for (i=0; i<size; i++) {
      if (count[i]<minload) minload = count[i];
      if (count[i]>maxload) maxload = count[i];
    }
    fprintf(stdout,
	    "#Load-TUW (min,exp,max): %d %d %d Imbalance %.2f\n",minload,ncut,maxload,100.0*((double)(maxload-minload)/ncut));

    if (!weak) {
      assert(m==n);
      d = (base_t*)malloc(m*sizeof(base_t));
      assert(d!=NULL);

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
      
      assert(memcmp(a,d,m*sizeof(base_t))==0);
      free(d);
    }
  } else {
    if (!weak) {
      assert(displ==NULL); assert(count==NULL); // just make sure
      d = NULL;

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
    }
  }
  if (c!=b) free(c); 
#endif

#ifdef WUTQSORT
  for (r=0; r<REPEAT; r++) partime[r] = 0.0;
  c = NULL;
  MPI_Barrier(MPI_COMM_WORLD);
  for (r=0,rr=0; r<REPEAT+WARMUP; r++) {
    memcpy(b,aa,np*sizeof(base_t));
    if (c!=b) free(c); // keep last
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    start = MPI_Wtime();
    mp = 0;
    srandom(seed);

    c = wut_qsort(b,np,&mp,MPI_COMM_WORLD);

    stop = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    if (r<WARMUP) continue;
    partime[rr++] = stop-start;
  }
  MPI_Allreduce(MPI_IN_PLACE,partime,REPEAT,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

  if (rank==0) {
    double paravg, parmin;

    paravg = 0.0; 
    parmin = partime[0];
    for (r=0; r<REPEAT; r++) {
      paravg += partime[r]; 
      if (partime[r]<parmin) parmin = partime[r];
    }
    paravg /= REPEAT;

    fprintf(stdout,
	    "P %d NP %d Seq %.2f WUT %.2f (avg) %.2f (min)\n",
	    size,np,seqtime*MICRO,paravg*MICRO,parmin*MICRO);
    fprintf(stdout,"#SU-WUT %.2f (avg) %.2f (min)\n",
	    seqtime/paravg,seqtime/parmin);
  }

  /* sorted? */
  sorted(c,mp,MPI_COMM_WORLD);
  
  /* compare against globally sorted */

  MPI_Gather(&mp,1,MPI_INT,count,1,MPI_INT,0,MPI_COMM_WORLD);
  if (rank==0) {
    displ[0] = 0;
    for (i=1; i<size; i++) displ[i] = displ[i-1]+count[i-1];
    m = displ[size-1]+count[size-1];

    int minload, maxload;
    minload = n; maxload = 0;
    for (i=0; i<size; i++) {
      if (count[i]<minload) minload = count[i];
      if (count[i]>maxload) maxload = count[i];
    }
    fprintf(stdout,
	    "#Load-WUT (min,exp,max): %d %d %d Imbalance %.2f\n",minload,ncut,maxload,100.0*((double)(maxload-minload)/ncut));

    if (!weak) {
      assert(m==n);
      d = (base_t*)malloc(m*sizeof(base_t));
      assert(d!=NULL);

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
      
      assert(memcmp(a,d,m*sizeof(base_t))==0);
      free(d);
    }
  } else {
    if (!weak) {
      assert(displ==NULL); assert(count==NULL); // just make sure
      d = NULL;

      MPI_Gatherv(c,mp,BASE_T,d,count,displ,BASE_T,0,MPI_COMM_WORLD);
    }
  }
  if (c!=b) free(c); 
#endif

  free(aa);
  free(b);

  if (rank==0) {
    free(a); 

    free(count);
    free(displ);
  }

  MPI_Op_free(&TUW_CHOOSE);
  MPI_Op_free(&TUW_SEGREP);

  MPI_Finalize();

  return 0;
}
