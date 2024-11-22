#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>


#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Custom reduction function
void custom_sum_function(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    int *in = (int *)invec;
    int *inout = (int *)inoutvec;
    for (int i = 0; i < *len; i++) {
        inout[i] += in[i];
    }
}

int main(int argc, char *argv[]) {
    int rank, size, ierr;
    int *send_buffer, *recv_buffer;
    MPI_Aint base_address, recv_disp;
    MPI_Aint displacements[1];
    int block_lengths[1];
    MPI_Datatype recv_type;
    MPI_Op custom_op;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Allocate and initialize buffers
    send_buffer = (int *)malloc(size * sizeof(int));
    recv_buffer = (int *)malloc(size * sizeof(int));

    for (int i = 0; i < size; i++) {
        send_buffer[i] = rank + 1; // Each process has unique values
        recv_buffer[i] = 0;        // Initialize receive buffer
    }

    // Get base address of recv_buffer
    MPI_Get_address(recv_buffer, &base_address);
    MPI_Get_address(&recv_buffer[0], &recv_disp);

    // Calculate displacement for recv_buffer relative to base address
    displacements[0] = recv_disp - base_address;

    // Define block length
    block_lengths[0] = size; // Number of integers in recv_buffer

    // Create a custom MPI datatype for the receive buffer
    MPI_Type_create_hindexed(1, block_lengths, displacements, MPI_INT, &recv_type);
    MPI_Type_commit(&recv_type);

    // Create a custom reduction operation
    MPI_Op_create(custom_sum_function, /*commute=*/1, &custom_op);

    // Perform an Allreduce operation using MPI_BOTTOM only for recv_buffer
    MPI_Allreduce(send_buffer, MPI_BOTTOM, 1, recv_type, custom_op, MPI_COMM_WORLD);

    // Print results
    printf("Process %d - send_buffer: ", rank);
    for (int i = 0; i < size; i++) {
        printf("%d ", send_buffer[i]);
    }
    printf("\n");

    printf("Process %d - recv_buffer after MPI_Allreduce: ", rank);
    for (int i = 0; i < size; i++) {
        printf("%d ", recv_buffer[i]);
    }
    printf("\n");

    // Clean up
    MPI_Type_free(&recv_type);
    MPI_Op_free(&custom_op);
    free(send_buffer);
    free(recv_buffer);
    MPI_Finalize();

    return 0;
}

