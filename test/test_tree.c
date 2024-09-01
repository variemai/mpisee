#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define N 56
typedef struct _tnode{
    char name[64];
    struct _tnode *parent;
    struct _tnode *child;
    struct _tnode *lsib;
    struct _tnode *rsib;
    struct _tnode **children;
    int nchildren;
    int bytes;
    int msgs;
}node;

int main(void){
    char names[N][16] = {"WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.0","WORLD","s_1.1","s_2.1","s_3.0","WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.1","WORLD","s_1.1","s_2.1","s_3.1","WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.0","WORLD","s_1.1","s_2.1","s_3.0","WORLD","s_1.0","s_2.0","WORLD","s_1.1","s_2.0","WORLD","s_1.0","s_2.1","s_3.1","WORLD","s_1.1","s_2.1","s_3.1"};
    char parents[N][16] = {"NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1","NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1","NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1","NULL","WORLD","s_1.0","NULL","WORLD","s_1.1","NULL","WORLD","s_1.0","s_2.1","NULL","WORLD","s_1.1","s_2.1"};
    uint64_t bytes[N] = {32,64,8,0,16,4,32,0,0,4,0,0,4,0,32,0,0,0,16,0,32,0,0,4,0,0,0,8,32,0,0,0,16,4,32,0,0,0,0,0,4,0,32,0,0,0,16,0,32,0,0,0,0,0,0,0};
    uint32_t msgs[N]= {1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0};
    int i,j,k;
    node *head = NULL;
    head = (node*) malloc (sizeof(node));
    strcpy(head->name, "WORLD");
    head->parent = NULL;
    head->bytes = 0;
    head->msgs = 0;
    head->lsib = NULL;
    head->rsib = NULL;
    head->child = NULL;
    node *new;
    for ( i = 0; i < N; i++){
        if ( strcmp(names[i],  "WORLD") == 0){
            head->bytes += bytes[i];
            head->msgs += msgs[i];
        }
    }
    k = 0;
    head->nchildren = 0;
    head->children = (node**) malloc ( sizeof(node*) *8 );
    for ( i =0; i<8; i++ ){
        head->children[i] = NULL;
    }
    node *ptr, *tmp;
    int found = 0;
    for ( i =1; i<N; i++ ){
        if ( strcmp(names[i],  "WORLD") != 0){
            new = (node*) malloc (sizeof(node));
            strcpy(new->name, names[i]);
            new->parent = NULL;
            ptr=head;
            /* find parent IT MUST ALREADY EXIST */
            found = 0;
            while ( ptr != NULL ){
                if ( strcmp(ptr->name, parents[i]) == 0 ){
                    found = 1;
                    break;
                }
                for ( i =0; i<ptr->nchildren; i++ ){
                    if ( strcmp(ptr->name, parents[i]) == 0 ){
                        found = 1;
                        break;
                    }
                }
                if (!found)

            }
            while ( ptr != NULL ){
                if ( strcmp(ptr->name, parents[i]) == 0 )
                    break;
                ptr=ptr->child;
            }
            if ( ptr == NULL  ){
                printf("should ptr be null?\n");
            }
            else{
                new->parent = ptr;
                printf("Parent of %s is %s\n",new->name,new->parent->name);
                ptr->children[ptr->nchildren] = new;
                ptr->nchildren = ptr->nchildren + 1;
                if ( ptr->child == NULL )
                    ptr->child = new;
                /* else{ */
                /*     ptr=ptr->child; */
                /*     while( ptr->rsib != NULL ) */
                /*         ptr = ptr->rsib; */
                /*     ptr->rsib = new; */
                /*     new->lsib = ptr->rsib; */
                /* } */
            }
            new->bytes = 0;
            new->msgs = 0;
        }
    }
    ptr = head;
    /* while(ptr != NULL){ */
    /*     printf("%s\n",ptr->name); */
    /*     ptr = ptr->child; */
    /* } */

    /* if ( strcmp(names[i],"WORLD") == 0 && head == NULL){ */
    /* } */
    return 0;
}
