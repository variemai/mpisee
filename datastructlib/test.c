#include "list.h"
#include <stdio.h>

void applyPrint(void **ptr,void *cl){
    printf("%s\n",(char*)*ptr);
}

int main(void){
    List_T lista;
    lista=List_list("arxodoa", "malakies", "kala", NULL);
    printf("%d\n",List_length(lista));
    List_map(lista,applyPrint,NULL);
    return 0;
}
