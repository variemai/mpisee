#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "except.h"
#include <stdint.h>
#define T Except_T

Except_Frame *Except_stack=NULL;

void Except_raise(const T *e,const char *file,int32_t line){
    Except_Frame *p = Except_stack;
    assert(e);
    if(p==NULL){
        fprintf(stderr,"Uncaught Exception");
        if(e->reason)fprintf(stderr, " %s",e->reason );
        else fprintf(stderr, " at 0x%p", (void*)e);
        if(file && line > 0) fprintf(stderr, " raised at %s:%d\n",file,line );
        fprintf(stderr, " Aborting...\n");
        fflush(stderr);
        abort();
    }
    p->exception = e;
    p->file = file;
    p->line = line;
    Except_stack=Except_stack->prev;
    longjmp(p->env,Except_raised);
}
