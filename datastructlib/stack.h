#ifndef _STACK_INCLUDED
#define _STACK_INCLUDED
#include <stdint.h>
typedef struct Stack_T *Stack_T;

extern Stack_T Stack_new(void);
extern int32_t Stack_empty(Stack_T stk);
extern void Stack_push(Stack_T stk,void *x);
extern void *Stack_pop(Stack_T stk);
extern void Stack_free(Stack_T *stk);

#endif
