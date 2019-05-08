#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>

#include <mpi.h>

#define MAXN 2000 /* Max value of N */
int N; /* Matrix size */
#define DIVFACTOR 32768.0f

#define SOURCE 0

/* My process rank           */
int my_rank;
/* The number of processes   */
int numProcess;

/* Matrixes given by a pointer */
float *A, *B, *X;



/* returns a seed for srand based on the time */
unsigned int time_seed() {
  struct timeval t;
  struct timezone tzdummy;
  gettimeofday(&t, &tzdummy);
  return (unsigned int)(t.tv_usec);
}


/* Set the program parameters from the command-line arguments */
void parameters(int argc, char **argv) {
  int seed = 0;  /* Random seed */
  char uid[32]; /*User name */

  /* Read command-line arguments */
  srand(time_seed());  /* Randomize */

  if (argc == 3) {
    seed = atoi(argv[2]);
    srand(seed);
    if ( my_rank == SOURCE ) printf("Random seed = %i\n", seed);
  }
  if (argc >= 2) {
    N = atoi(argv[1]);
    if (N < 1 || N > MAXN) {
      if ( my_rank == SOURCE ) printf("N = %i is out of range.\n", N);
      exit(0);
    }
  }
  else {
    if ( my_rank == SOURCE ) printf("Usage: %s <matrix_dimension> [random seed]\n",
           argv[0]);
    exit(0);
  }

  /* Print parameters */
  if ( my_rank == SOURCE )  printf("\nMatrix dimension N = %i.\n", N);
}


/* Allocates memory for A, B and X */
void allocate_memory() {
    A = (float*)malloc( N*N*sizeof(float) );
    B = (float*)malloc( N*sizeof(float) );
    X = (float*)malloc( N*sizeof(float) );
}

/* Free allocated memory for arrays */
void free_memory() {
    free(A);
    free(B);
    free(X);
}

/* Initialize A and B (and X to 0.0s) */
void initialize_inputs() {

  int row, col;

  for (row = 0; row < N; row++) {
    for (col = 0; col < N; col++) {
      A[col+N*row] = (float)rand() / DIVFACTOR;
    }
    B[row] = (float)rand() / DIVFACTOR;
    X[row] = 0.0;
  }

}


int main(int argc, char **argv) {

    /* Prototype functions*/
    void gauss();

    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* Every process reads the parameters to prepare dimension */
    parameters(argc, argv);

    /* Every process must allocate memory for the arrays */
    allocate_memory();

    if ( my_rank == SOURCE ) {
        /* Initialize A and B */
        initialize_inputs();
    }

    gauss();

    /* The barrier prevents any process to reach the finalize before the others have finished their communications */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Free memory used for the arrays that we allocated previously */
    free_memory();

    MPI_Finalize();
}


/* Includes both algorithms */
void gauss() {

    void gaussElimination();
    void backSubstitution();

    /* Times */
    double t1, t2;

    /* Barrier to sync all processes before starting the algorithms */
    MPI_Barrier(MPI_COMM_WORLD);

    /* Initial time */
    if ( my_rank == SOURCE )
        t1 = MPI_Wtime();

    /* Gauss Elimination is performed using MPI */
    gaussElimination();

    /* Back Substitution is performed sequentially */
    if ( my_rank == SOURCE ) {
        backSubstitution();

        /* Finish time */
        t2 = MPI_Wtime();

        printf("\nElapsed time: %f miliseconds\n", (t2-t1) * 1000 );
    }

}

/* Guassian Elimination algorithm using MPI */
void gaussElimination() {

    MPI_Status status;
    MPI_Request request;
    int row, col, i, norm;
    float multiplier;

    /* Array with the row size and number of rows that each processor will handle */
    int * first_row_A_array = (int*) malloc ( p * sizeof(int) );
    int * n_of_rows_A_array = (int*) malloc ( p * sizeof(int) );
    int * first_row_B_array = (int*) malloc ( p * sizeof(int) );
    int * n_of_rows_B_array = (int*) malloc ( p * sizeof(int) );
    for ( i = 0; i < p; i++ ) {
        first_row_A_array[i] = 0;
        n_of_rows_A_array[i] = 0;
        first_row_B_array[i] = 0;
        n_of_rows_B_array[i] = 0;
    }

    /* Main loop. After every iteration, a new column will have all 0 values down the [norm] index */
    for (norm = 0; norm < N-1; norm++) {

        /* Broadcast the A[norm] row and B[norm], important values of this iteration */
        MPI_Bcast( &A[ N*norm ], N, MPI_FLOAT, SOURCE, MPI_COMM_WORLD );
        MPI_Bcast( &B[norm], 1, MPI_FLOAT, SOURCE, MPI_COMM_WORLD );

        /* subset of rows of this iteration */
        int subset = N - 1 - norm;
        /* number that indicates the step as a float */
        float step = ((float)subset ) / (numProcess);
        /* First and last rows that this process will work into for this iteration */
        int first_row = norm + 1 + ceil( step * (my_rank) );
        int last_row = norm + 1 + floor( step * (my_rank+1) );
        if ( last_row >= N ) last_row = N-1;
        int number_of_rows = last_row - first_row +1;

        if ( my_rank == SOURCE ) {

            for ( i = 1; i < p; i++ ) {

                /* We send to each process the amount of data that they are going to handle */
                int first_row_rmte = norm + 1 + ceil( step * (i) );
                int last_row_rmte = norm + 1 + floor( step * (i+1) );
                if( last_row_rmte >= N ) last_row_rmte = N -1;
                int number_of_rows_rmte = last_row_rmte - first_row_rmte +1;

                if ( number_of_rows_rmte < 0 ) number_of_rows_rmte = 0;
                if ( first_row_rmte >= N ) { number_of_rows_rmte = 0; first_row_rmte = N-1; };

                first_row_A_array[i] = first_row_rmte * N;
                first_row_B_array[i] = first_row_rmte;
                n_of_rows_A_array[i] = number_of_rows_rmte * N;
                n_of_rows_B_array[i] = number_of_rows_rmte ;

            }

        }

        MPI_Scatterv(
            &A[0],              // send buffer
            n_of_rows_A_array,  // array with number of elements in each chunk
            first_row_A_array,  // array with pointers to initial element of each chunk
            MPI_FLOAT,          // type of elements to send
            &A[first_row * N],  // receive buffer
            N * number_of_rows, // number of elements to receive
            MPI_FLOAT,          // type of elements to receive
            SOURCE,             // who sends
            MPI_COMM_WORLD
        );
        MPI_Scatterv(
            &B[0],
            n_of_rows_B_array,
            first_row_B_array,
            MPI_FLOAT,
            &B[first_row],
            number_of_rows,
            MPI_FLOAT,
            SOURCE,
            MPI_COMM_WORLD
        );

        if ( number_of_rows > 0  && first_row < N) {

            for (row = first_row; row <= last_row; row++) {

                multiplier = A[N*row + norm] / A[norm + N*norm];
                for (col = norm; col < N; col++) {
                    A[col+N*row] -= A[N*norm + col] * multiplier;
                }

                B[row] -= B[norm] * multiplier;
            }
        }

        if ( my_rank != SOURCE ) {
            if ( number_of_rows > 0  && first_row < N) {
                MPI_Isend( &A[first_row * N], N * number_of_rows, MPI_FLOAT, SOURCE,0, MPI_COMM_WORLD, &request);
                MPI_Isend( &B[first_row],         number_of_rows, MPI_FLOAT, SOURCE,0, MPI_COMM_WORLD, &request);
            }
        }

        else {

            for ( i = 1; i < p; i++ ) {

                if( n_of_rows_B_array[i] < 1  || first_row_B_array[i] >= N) continue;

                MPI_Recv( &A[ first_row_A_array[i] ], n_of_rows_A_array[i] , MPI_FLOAT, i,0, MPI_COMM_WORLD, &status );
                MPI_Recv( &B[ first_row_B_array[i] ], n_of_rows_B_array[i] , MPI_FLOAT, i,0, MPI_COMM_WORLD, &status );
            }


        }

    }
}


/* Back substitution sequential algorithm */
void backSubstitution () {

    int norm, row, col;  /* Normalization row, and zeroing
      * element row and col */

    /* Back substitution */
    for (row = N - 1; row >= 0; row--) {
        X[row] = B[row];

        for (col = N-1; col > row; col--) {
            X[row] -= A[N*row+col] * X[col];
        }

        X[row] /= A[N*row + row];
    }
}
