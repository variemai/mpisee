#include <assert.h>
#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    int buffer[1];                  // Buffer for data to broadcast
    MPI_Aint addr;                  // Base address for buffer
    int ones[1] = {1};              // Block length for derived datatype
    MPI_Datatype tint[1] = {MPI_INT}; // Element type (MPI_INT)
    MPI_Datatype bcast_type;        // Derived datatype

    MPI_Init(&argc, &argv);

    // Get the base address of the buffer
    MPI_Get_address(buffer, &addr);

    // Create a derived datatype for the buffer using MPI_BOTTOM
    MPI_Type_create_struct(1, ones, &addr, tint, &bcast_type);
    MPI_Type_commit(&bcast_type);

    // Initialize the buffer differently on rank 0
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        buffer[0] = 42; // Root initializes buffer
    } else {
        buffer[0] = -1; // Other ranks initialize buffer to a dummy value
    }

    // Broadcast the data using MPI_BOTTOM
    MPI_Bcast(MPI_BOTTOM, 1, bcast_type, 0, MPI_COMM_WORLD);

    // Verify the result
    assert(buffer[0] == 42);
    printf("Rank %d: buffer[0] = %d\n", rank, buffer[0]);

    // Clean up
    MPI_Type_free(&bcast_type);
    MPI_Finalize();

    return 0;
}
