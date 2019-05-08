#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

#ifndef MPI_WRAPPERS_H
#define	MPI_WRAPPERS_H

#ifdef	__cplusplus
extern "C" {
#endif

    void wait_all_wrapper(MPI_Request* requests, int size);
    void wait_wrapper(MPI_Request* request);

#ifdef	__cplusplus
}
#endif

#endif	/* MPI_WRAPPERS_H */

