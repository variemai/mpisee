#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char *argv[]) {
    int world_rank, world_size;
    int local_rank, local_size;
    MPI_Comm comm_local;
    int color;
    MPI_Request requests[2]; // 0: world request, 1: local request
    MPI_Status statuses[2];

    int send_buf_world, recv_buf_world;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Check if we have exactly 2 processes
    if (world_size != 2) {
        if (world_rank == 0) {
            fprintf(stderr, "This program requires exactly 2 MPI processes.\n");
        }
        MPI_Finalize();
        return 1;
    }

    color = world_rank % 2;

    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &comm_local);

    // Get rank and size within the new local communicator
    MPI_Comm_rank(comm_local, &local_rank);
    MPI_Comm_size(comm_local, &local_size);

    // --- Communication Phase ---

    send_buf_world = world_rank; // Simple data
    MPI_Iallreduce(&send_buf_world, &recv_buf_world, 1, MPI_INT, MPI_SUM,
                   MPI_COMM_WORLD, &requests[0]);
    MPI_Ibcast(&send_buf_world, 1, MPI_INT, 0, comm_local, &requests[1]);


    // --- Wait Phase ---
    // Now wait for BOTH requests to complete.
    // This is where MPisee faces the challenge: requests[0] is associated
    // with MPI_COMM_WORLD, requests[1] with comm_local.
    printf("World Rank %d: Calling MPI_Waitall on 2 requests...\n", world_rank);
    MPI_Waitall(2, requests, statuses);
    printf("World Rank %d: MPI_Waitall completed.\n", world_rank);

    // --- Cleanup ---
     MPI_Comm_free(&comm_local);
    MPI_Finalize();

    return 0;
}
