/**
 * @author RookieHPC
 * @brief Original source code at https://www.rookiehpc.com/mpi/docs/mpi_cart_sub.php
 **/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>
#include <string.h>

/**
 * @brief Illustrates how to partition a cartesian topology created with
 * MPI_Cart_create.
 * @details This program is meant to be run with 6 MPI processes. It consists in
 * partitioning a 2 x 3 2D cartesian topology along its first dimension to
 * obtain 2 1D subgrids of 3 MPI processes each.
 * The initial 2D cartesian topology can be visualised as follows:
 *
 *                  Second dimension
 * ------------------------------------------------->
 *
 * +---------------+---------------+---------------+   |
 * | MPI process 0 | MPI process 1 | MPI process 2 |   |
 * |     (0,0)     |     (0,1)     |     (0,2)     |   |
 * +---------------+---------------+---------------+   |  First dimension
 * | MPI process 3 | MPI process 4 | MPI process 5 |   |
 * |     (1,0)     |     (1,1)     |     (1,2)     |   |
 * +---------------+---------------+---------------+   v
 *
 * And the final subgrids can be visualised as follows:
 *
 * +---------------+---------------+---------------+  \
 * | MPI process 0 | MPI process 1 | MPI process 2 |   } First subgrid
 * +---------------+---------------+---------------+  /
 *
 * +---------------+---------------+---------------+  \
 * | MPI process 3 | MPI process 4 | MPI process 5 |   } Second subgrid
 * +---------------+---------------+---------------+  /
 *
 **/
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    // Size of the default communicator
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(size != 24)
    {
        printf("This application is meant to be run with 24 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Ask MPI to decompose our processes in a 2D cartesian grid for us
    int dims[3] = {3, 4, 2};

    // Make both dimensions non-periodic
    int periods[3] = {false, false ,false};

    // Let MPI assign arbitrary ranks if it deems it necessary
    int reorder = true;

    // Create a communicator given the 2D torus topology.
    MPI_Comm cartesian_communicator;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, reorder, &cartesian_communicator);

    // My rank in the new communicator
    int my_rank;
    MPI_Comm_rank(cartesian_communicator, &my_rank);

    // Get my coordinates in the new communicator
    int my_coords[3];
    MPI_Cart_coords(cartesian_communicator, my_rank, 3, my_coords);

    // Print my location in the 2D cartesian topology.
    printf("[MPI process %d] I am located at (%d, %d, %d) in the initial 2D cartesian topology.\n", my_rank, my_coords[0], my_coords[1], my_coords[2]);

    // Partition the 2D cartesian topology along the first dimension, by preserving the second dimension
    int remain_dims[3] = {0, 1, 0};
    MPI_Comm subgrid_communicator;
    MPI_Cart_sub(cartesian_communicator, remain_dims, &subgrid_communicator);
    int *subgrid_ranks, ssg_size;
    MPI_Comm_size(subgrid_communicator, &ssg_size);
    subgrid_ranks = (int*) malloc (sizeof(int)*ssg_size);

    // Get the ranks of all MPI processes in my subgrid and print it
    MPI_Allgather(&my_rank, 1, MPI_INT, subgrid_ranks, 1, MPI_INT, subgrid_communicator);
    /* printf("[MPI process %d] I am in the 1D subgrid that contains MPI processes %d, %d %d and %d.\n", my_rank, subgrid_ranks[0], subgrid_ranks[1], subgrid_ranks[2],subgrid_ranks[3]); */
    int i,j,newcomms,maxdims,subdims[6],subperiods[6],subcoords[6],realdims;
    int *all_ranks,*uniq;
    newcomms = 1;
    maxdims = 6;
    realdims = 0;
    memset(subdims, 0, sizeof(int)*6);
    MPI_Cart_get(cartesian_communicator, maxdims, subdims, subperiods, subcoords);
    if ( my_rank == 0 ){
            for ( i =0; i < maxdims; i++ ){
                if ( subdims[i] != 0 )
                    realdims++;
            }
            printf("Realdims %d %d %d %d\n",realdims,subdims[0],subdims[1],subdims[2]);
            //printf("MAXDIMS  = %d\n",maxdims);
    }
    for ( i =0; i<maxdims; i++ ){
        if ( !remain_dims[i]  )
            newcomms = newcomms * subdims[i];
    }
    if ( my_rank == 0 )
        printf("Newcoms = %d\n",newcomms);
    /* all_ranks = (int*) malloc (size*sizeof(int)); */
    /* uniq = (int*) calloc (newcomms,sizeof(int)); */
    /* MPI_Gather(&subgrid_ranks[0], 1, MPI_INT, all_ranks,1, MPI_INT, 0, MPI_COMM_WORLD); */
    /* if ( my_rank == 0 ){ */
    /*     for ( i =0 ; i<size; i++ ){ */
    /*         printf("Recvd rank = %d\n",all_ranks[i]); */
    /*     } */
    /* } */
    /* int found; */
    /* int k = 0; */
    /* if ( my_rank == 0 ){ */
    /* for (i = 0; i < size; i++) { */
    /*     found = 0; */
    /*     for ( j =0; j<newcomms; j++ ){ */
    /*         if ( all_ranks[i] == uniq[j] ){ */
    /*             found = 1; */
    /*             break; */
    /*         } */
    /*     } */
    /*     if ( !found  ){ */
    /*         uniq[k] = all_ranks[i]; */
    /*         k++; */
    /*     } */
    /* } */
    /*     for ( i =0; i<newcomms; i++ ){ */
    /*         printf("%d\n",uniq[i]); */
    /*     } */
    /* } */

    MPI_Finalize();

    return EXIT_SUCCESS;
}
