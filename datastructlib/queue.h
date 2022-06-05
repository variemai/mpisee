#ifndef _QUEUE_INCLUDED
#define _QUEUE_INCLUDED
#include <stdint.h>
typedef struct Queue_T *Queue_T;

extern Queue_T Queue_new(void);
extern int32_t Queue_empty(Queue_T que);
extern void Queue_enqueue(Queue_T que,void *x);
extern void *Queue_dequeue(Queue_T que);
extern void Queue_free(Queue_T *que);

#endif
