
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <unistd.h>

#include "gauss_jordan.h"

gauss_jordan gj;

int main(int argc, char** argv) {
    srand(time(NULL));
    if (argc != 3) {
        fprintf(stderr, "Usage: %s [dimension (not augmented)] [use k-group distribution]\n",
                argv[0]);
        exit(EXIT_SUCCESS);
    }
    int dimension = atoi(argv[1]);
    if (dimension <= 0 || (dimension % 2 != 0)) {
        fprintf(stderr, "[dimension] must be a positive even number\n");
        exit(EXIT_SUCCESS);
    }
    int groupsDistribution = atoi(argv[2]);
    if (groupsDistribution < 0 || groupsDistribution > 1) {
        fprintf(stderr, "[use k-group distribution must be either 0 or 1]\n");
        exit(EXIT_SUCCESS);
    }
    // Create augmented random array
    float** augmented_m = create_augmented_matrix(dimension);
    int augmented_n = dimension + 1;


    column_t* my_cols = init(augmented_n, augmented_m, groupsDistribution);
    MPI_Barrier(MPI_COMM_WORLD);
    gj_kgi_main_loop(my_cols);
    if (gj.my_rank == 0) {

    }

    MPI_Barrier(MPI_COMM_WORLD);
    destroy(&my_cols);
    MPI_Finalize();
    return (EXIT_SUCCESS);
}

