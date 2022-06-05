#ifndef TABLE_INCLUDED
#define TABLE_INCLUDED

#include <assert.h>
#include <stddef.h>
#include "mem.h"
#include <stdint.h>
#include "atom.h"

#define T Table_T
typedef struct T *T;

extern T Table_new(int32_t hint32_t, int32_t cmp(const void *x,const void *y),
                   unsigned hash(const void *key));
extern void Table_free(T* table);
extern int Table_length(T table);
extern void* Table_put(T table, const void *key, void *value);
extern void* Table_get(T table, const void* key);
extern void* Table_remove(T table, const void* key);
extern void Table_map(T table, void apply(const void *key, void **value,
                                          void *cl),void *cl);
extern void** Table_toArray(T table, void* end);

#endif
