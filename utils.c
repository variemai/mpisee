#include "utils.h"
#include "symbols.h"
#define MAX_ARG_STRING_SIZE 4096

const char prim_names[][NUM_OF_PRIMS]={
"Send",
"Recv",
"Isend",
"Irecv",
"Ssend",
"Issend",
"Probe",
"Iprobe",
"Waitall",
"Wait",
"Waitany",
"Test",
"Testany",
"Testall",
"Sendrecv",
"Bcast",
"Barrier",
"Allreduce",
"Allgather",
"Allgatherv",
"Alltoall",
"Alltoallv",
"Alltoallw",
"Reduce",
"Gather",
"Gatherv",
"Scan",
"Exscan",
"Scatter",
"Scatterv",
"Reduce_scatter",
"Iallreduce",
"Ibcast",
"Ialltoall",
"Iscatter",
"Ibarrier"
};

char *appname = NULL;
prof_attrs **comm_table = NULL;

void mcpt_abort (char *fmt, ...){
  va_list args;
  va_start (args, fmt);
  fprintf (stderr, "\n\n: MPICP ABORTING ");
  vfprintf (stderr, fmt, args);
  va_end (args);
  fflush (stderr);
  PMPI_Abort(MPI_COMM_WORLD, -1);
}

char *get_appname (){
  int pid, exelen, insize = 256;
  char *inbuf = NULL, file[256];

  pid = getpid ();
  snprintf (file, 256, "/proc/%d/exe", pid);
  inbuf = malloc (insize);
  if (inbuf == NULL){
      mcpt_abort ("unable to allocate space for full executable path.\n");
  }

  exelen = readlink (file, inbuf, 256);
  if (exelen == -1){
      if (errno != ENOENT){
          while (exelen == -1 && errno == ENAMETOOLONG){
              insize += 256;
              inbuf = realloc (inbuf, insize);
              exelen = readlink (file, inbuf, insize);
          }
          inbuf[exelen] = '\0';
          return inbuf;
      }
      else
        free (inbuf);
  }
  else{
      inbuf[exelen] = '\0';
      return inbuf;
  }
  return NULL;
}

void
getProcCmdLine (int *ac, char **av)
{
  int i = 0, pid;
  char *inbuf, file[256];
  FILE *infile;
  char *arg_ptr;

  *ac = 0;
  *av = NULL;

  pid = getpid ();
  snprintf (file, 256, "/proc/%d/cmdline", pid);
  infile = fopen (file, "r");

  if (infile != NULL){
    while (!feof (infile)){
      inbuf = malloc (MAX_ARG_STRING_SIZE);
      if (fread (inbuf, 1, MAX_ARG_STRING_SIZE, infile) > 0){
        arg_ptr = inbuf;
        while (*arg_ptr != '\0'){
          av[i] = strdup (arg_ptr);
          arg_ptr += strlen (av[i]) + 1;
          i++;
        }
      }
    }
    *ac = i;
    if(inbuf != NULL ){
      free (inbuf);
    }
    fclose (infile);
  }
  else{
    mcpt_abort("Error opening file %s FILE:LINE = %d",file,__FILE__,__LINE__);
  }
}
