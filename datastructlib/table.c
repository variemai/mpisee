#include "table.h"
#include <limits.h>
#include "mem.h"
#include <stdio.h>

struct T{
    int size;
    int length;
    unsigned timestamp;
    int (*cmp)(const void *x, const void *y);
    unsigned (*hash)(const void *key);
    struct binding {
        struct binding *link;
        const void *key;
        void *value;
    }**buckets;
};

static int _cmpatom(const void *x, const void *y) {
    return x != y;
}

static unsigned _hashatom(const void *key) {
    return (unsigned long)key>>2;
}

T Table_new(int32_t hint, int32_t cmp(const void *x,const void *y),
          unsigned hash(const void *key)){
    T table;
    int i;
    static int primes [] = {251, 509, 1021, 2053, 4093, 8191, 16381, 32771,
    65521,127913, INT_MAX};
    assert( hint >= 0);
    for ( i = 1; primes[i] < hint; i++ )
        ;
    table = ALLOC(sizeof (*table) + primes[i-1]*sizeof(table->buckets[0]));
    table->size = primes[i-1];
    table->cmp = cmp ? cmp : _cmpatom;
    table->hash = hash ? hash : _hashatom;
    table->buckets = (struct binding**)(table+1);
    for (i =0; i<table->size; i++)
        table->buckets[i] = NULL;
    table->length = 0;
    table->timestamp = 0;
    return table;
}

void* Table_put(T table, const void* key, void* value){
    int i;
    struct binding *p;
    void *prev;
    assert(table);
    assert(key);
    i= (*table->hash)(key)%table->size;
    for ( p = table->buckets[i]; p; p = p->link ){
        if ( (*table->cmp)(key, p->key) == 0 )
            break;
    }
    if ( p == NULL ){
        NEW(p);
        p->key = key;
        p->link = table->buckets[i];
        table->buckets[i] = p;
        table->length++;
        prev = NULL;
    }
    else{
        prev = p->value;
    }
    p->value = value;
    table->timestamp++;
    return prev;
}


/*****************************************************************************/
/* Table_get finds a binding by hashing its key, taking it modulo the number */
/* of elements in buckets, and searching the list for a key equal to key.    */
/* It calls the table’s hash and cmp functions. This for loop terminates     */
/* when it finds the key, and it thus leaves p pointing to the binding of    */
/* interest. Otherwise, p ends up null                                       */
/*****************************************************************************/
void* Table_get(T table, const void* key){
    int i;
    struct binding *p = NULL;
    assert(table);
    if ( !key )
        return NULL;
    /* assert(key); */
    i = (*table->hash)(key)%table->size;
    for ( p = table->buckets[i]; p; p = p->link ){
        if ( ( *table->cmp )(key,p->key) == 0 ){
            break;
        }
    }
    return p ? p->value : NULL;
}

void Table_map(T table, void apply(const void *key, void **value, void *cl),
       void *cl) {
    int i;
    unsigned stamp;
    struct binding *p;

    assert(table);
    assert(apply);
    stamp = table->timestamp;
    for (i = 0; i < table->size; i++){
        for (p = table->buckets[i]; p; p = p->link) {
            apply(p->key, &p->value, cl);
            assert(table->timestamp == stamp);
        }
    }
}

/****************************************************************************/
/* Table_toArray allocates an array to hold the key-value pairs followed by */
/* a terminating end pointer, and fills in the array by visiting each       */
/* binding in table. p->key must be cast from const void * to void *        */
/* because the array is not declared const. The order of the key-value      */
/* pairs in the array is arbitrary.                                         */
/****************************************************************************/
void** Table_toArray(T table, void* end){
    int i,j;
    void **array;
    struct binding *p;
    assert(table);
    j = 0;
    array = ALLOC((2*table->length + 1)*sizeof(*array));
    if ( array == NULL ){
        fprintf(stderr,"ALLOC %d\n",__LINE__);
    }
    for (i = 0; i < table->size; i++) {
        for (p = table->buckets[i]; p; p=p->link) {
            array[j++] = (void*)p->key;
            array[j++] = p->value;
        }
    }
    array[j] = end;
    return array;
}

void *Table_remove(T table, const void *key) {
    int i;
    struct binding **pp;

    assert(table);
    assert(key);
    /* if ( !key ) */
    /*     return NULL; */
    table->timestamp++;
    i = (*table->hash)(key)%table->size;
    for (pp = &table->buckets[i]; *pp; pp = &(*pp)->link)
        if ((*table->cmp)(key, (*pp)->key) == 0) {
            struct binding *p = *pp;
            void *value = p->value;
            *pp = p->link;
            FREE(p);
            table->length--;
            return value;
        }
    return NULL;
}


void Table_free(T *table){
    assert(table && *table);
    if ( (*table)-> length > 0 ){
        int i;
        struct binding *p, *q;
        for ( i =0; i< (*table)->size; i++ ){
            for ( p =(*table)->buckets[i]; p; p=q ){
                q = p->link;
                FREE(p);
            }
        }
    }
    FREE(*table);
}

extern int Table_length(T table){
    assert(table);
    return table->length;
}
