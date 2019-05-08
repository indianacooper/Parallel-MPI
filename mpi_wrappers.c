#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mpi_wrappers.h"
#include "gauss_jordan.h"

extern gauss_jordan gj;
/**
 * Waiting for all requests array of length size
 * @param requests The MPI_Requests array
 * @param size The length of the array
 */
void wait_all_wrapper(MPI_Request* requests, int size) {
    MPI_Status* statuses = malloc(size * sizeof(MPI_Status));
    int wait_all_retval = MPI_Waitall(size, requests, statuses);
    if (wait_all_retval != MPI_SUCCESS) {
        switch(wait_all_retval) {
            case (MPI_ERR_REQUEST):
                fprintf(stderr, "wait_all_wrapper: MPI_ERR_REQUEST\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            case (MPI_ERR_ARG):
                fprintf(stderr, "wait_all_wrapper: MPI_ERR_ARG\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            case (MPI_ERR_IN_STATUS):
                fprintf(stderr, "wait_all_wrapper: MPI_ERR_IN_STATUS\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            default:
                MPI_Finalize();
                exit(EXIT_FAILURE);
        }
    }
    free(statuses);
    statuses = NULL;
}

/**
 * Waiting the given request
 * @param request The MPI_Request 
 */
void wait_wrapper(MPI_Request* request) {
    MPI_Status mpi_status;
    int status = MPI_Wait(request, &mpi_status);
    if (status != MPI_SUCCESS) {
        switch(status) {
            case (MPI_ERR_REQUEST):
                fprintf(stderr, "wait_wrapper: MPI_ERR_REQUEST\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            case (MPI_ERR_ARG):
                fprintf(stderr, "wait_wrapper: MPI_ERR_ARG\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            case (MPI_ERR_IN_STATUS):
                fprintf(stderr, "wait_wrapper: MPI_ERR_IN_STATUS\n");
                MPI_Finalize();
                exit(EXIT_FAILURE);
            default:
                MPI_Finalize();
                exit(EXIT_FAILURE);
        }
    }
}
