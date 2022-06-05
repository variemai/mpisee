#include "queue.h"
#include <assert.h>
#include <stddef.h>
#include "mem.h"
#include <stdint.h>

#define T Queue_T

struct elem{
    void *x;
    struct elem *link;
};

struct T{
    int32_t count;
    struct elem* head;
    struct elem* tail;
};

T Queue_new(void){
    T que;
    NEW(que);
    que->count = 0;
    que->head = NULL;
    que->tail = NULL;
    return que;
}

int32_t Queue_empty(T que){
    assert(que);
    return que->count == 0;
}

void Queue_enqueue(T que, void *x){
    struct elem *t;
    assert(que);
    NEW(t);
    t->x=x;
    t->link=NULL;
    if(Queue_empty(que)){
        que->head = t;
        que->tail = t;
    }
    else{
        que->tail->link = t;
        que->tail = t;
    }
    que->count++;
}

void *Queue_dequeue(T que){
    void *x;
    struct elem *t;
    assert(que);
    assert(que->count>0);
    t=que->head;
    que->head=t->link;
    que->count--;
    x=t->x;
    FREE(t);
    return x;
}

void Queue_free(T* que){
    struct elem *t, *u;
    assert(que && *que);
    for(t=(*que)->head; t; t=u){
        u=t->link;
        FREE(t);
    }
    FREE(*que);
}
