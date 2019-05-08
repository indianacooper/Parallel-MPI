#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "gauss_jordan.h"
#include "column.h"
#include "mpi_wrappers.h"

extern gauss_jordan gj;

// Forward declarations
void wait_all_wrapper(MPI_Request* requests, int size);
int is_my_column(int col_idx);
int process_of_column(int col_idx);
int number_of_cols();
int map_global_to_local(int k);
int number_of_recipients();
MPI_Request* send_column(int k, int max_idx, column_t* my_cols);
column_t receive_column(int k);
void destroy_augmented_matrix(float** augmented_m);
void count_processes();

/**
 * Initializes the program, MPI framework, and sends the columns to the
 * appropriate MPI processes
 * @param augmented_n The length of the 2D array (B-vector inclusive)
 * @param augmented_m The 2D array of coefficients
 * @param groupDistribution
 * @return The array of columns of each process
 */
column_t* init(int augmented_n, float** augmented_m, int useGroupDistribution) {
    gj.dimension = augmented_n - 1;
    gj.use_group_distribution = (useGroupDistribution != 0) ? 1 : 0;
    // Init MPI
    MPI_Init(NULL, NULL);
    // Initialize global Gauss Jordan info struct
    MPI_Comm_size(MPI_COMM_WORLD, &(gj.proc_num));
    if (gj.dimension % gj.proc_num != 0) {
        fprintf(stderr, "Dimension given, not divisible by process_number\n");
        exit(EXIT_FAILURE);
    }
    gj.group_number = gj.dimension / gj.proc_num;

    MPI_Comm_rank(MPI_COMM_WORLD, &gj.my_rank);
    column_t* my_cols;
    // Master populates an all-columns array
    // then sends the appropriate columns to their "owners"
    // Others, receive their columns and add them to their column array
    if (gj.my_rank == 0) {
        // Convert array to columns
        float** columns = malloc(augmented_n * sizeof(float*));
        int j;
        for (j = 0; j < augmented_n; j++) {
            columns[j] = malloc((augmented_n - 1) * sizeof(float));
            int i;
            for (i = 0; i < augmented_n - 1; i++) {
                columns[j][i] = augmented_m[i][j];
            }
        }

        // Send to others and populate yours
        // (master handles the b-vector as well)
        my_cols = malloc((gj.group_number + 1) * sizeof(column_t));
        memset(my_cols, '\0', (gj.group_number + 1) * sizeof(column_t));
        MPI_Request* requests =
            malloc(gj.group_number * (gj.proc_num - 1) * sizeof(MPI_Request));
        memset(requests, '\0', gj.group_number * (gj.proc_num - 1) * sizeof(MPI_Request));
        int my_cols_idx = 0, requests_idx = 0;
        for (j = 0; j < augmented_n - 1; j++) {
            int others_rank = process_of_column(j);
            if (others_rank != gj.my_rank)
                MPI_Isend(columns[j], gj.dimension, MPI_FLOAT, others_rank, j,
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
            else
                my_cols[my_cols_idx++] = create_column(j, columns[j]);
        }
        // Add b-vector
        my_cols[my_cols_idx] = create_column(j, columns[gj.dimension]);

        // Wait on sent data
        wait_all_wrapper(requests, requests_idx);
        free(requests);
        requests = NULL;
        // Free all columns array
        for (j = 0; j < augmented_n; j++) {
            free(columns[j]);
            columns[j] = NULL;
        }
        free(columns);
        columns = NULL;
    }
    else {
        // Allocate space for your columns and receive them (synchronously)
        my_cols = malloc(gj.group_number * sizeof(column_t));
        memset(my_cols, '\0', gj.group_number * sizeof(column_t));
        float* tmp_col = malloc(gj.dimension * sizeof(float));
        memset(tmp_col, '\0', gj.dimension * sizeof(float));
        int j, my_cols_idx = 0;
        for (j = 0; j < gj.dimension; j++) {
            if (is_my_column(j)) {
                MPI_Request request;
                MPI_Irecv(tmp_col, gj.dimension, MPI_FLOAT, 0,
                    j, MPI_COMM_WORLD, &request);
                wait_wrapper(&request);
                my_cols[my_cols_idx++] = create_column(j, tmp_col);
                if (my_cols_idx == gj.group_number)
                    break;
            }
        }
        free(tmp_col);
        tmp_col = NULL;
    }
    // Allocate space for the dummy column, used for sending and receiving
    // a column (acts as a serialization array)
    gj.dummy_col = malloc((gj.dimension + 1) * sizeof(float));
    destroy_augmented_matrix(augmented_m);
    return my_cols;
}

void destroy(column_t** my_cols_ptr) {
    // Free my_cols array
    int i;
    int upper_dim = (gj.my_rank == 0) ? gj.group_number + 1 : gj.group_number;
    for (i = 0; i < upper_dim; i++) {
        delete_column((*my_cols_ptr)[i]);
        (*my_cols_ptr)[i] = NULL;
    }
    free(*my_cols_ptr);
    *my_cols_ptr = NULL;
    // Free dummy column
    free(gj.dummy_col);
    gj.dummy_col = NULL;
}

void gj_kgi_main_loop(column_t* my_cols) {
    int k;
    double start_time = MPI_Wtime();
    for (k = 0; k < gj.dimension; k++) {
        if (is_my_column(k)) {
            column_t col = my_cols[map_global_to_local(k)];
            int pivot_idx = pivot(col, k);
            MPI_Request* requests = send_column(k, pivot_idx, my_cols);
            int j, upper_bound = number_of_cols();
            for (j = map_global_to_local(k) + 1; j < upper_bound; j++) {
                modify(my_cols[j], col, pivot_idx, k);
            }
            wait_all_wrapper(requests, gj.recipients_number);
            free(requests);
            requests = NULL;
        }
        else {
            // Receive column and extract pivot index
            if ((gj.my_rank == 0) ||
                    ( !(gj.use_group_distribution && (gj.my_rank < process_of_column(k))) )) {
                if (!gj.use_group_distribution && (k > gj.dimension - gj.proc_num) &&
                    (process_of_column(k) > gj.my_rank))
                    break;
                column_t col = receive_column(k);
                int pivot_idx = (int) gj.dummy_col[0];

                if (gj.use_group_distribution) {
                    if (gj.my_rank == 0)
                        modify(my_cols[gj.group_number], col, pivot_idx, k);
                    else {
                        int j;
                        for (j = 0; j < gj.group_number; j++)
                            modify(my_cols[j], col, pivot_idx, k);
                    }
                }
                else {
                    int local_idx = map_global_to_local(k);
                    // Not the latest column, start modifying from the local_idx
                    int j, start_idx = (process_of_column(k) != gj.proc_num - 1) ?
                            local_idx : local_idx + 1;
                    int upper_bound = (gj.my_rank == 0) ?
                        gj.group_number + 1 : gj.group_number;
                    for (j = start_idx; j < upper_bound; j++)
                        modify(my_cols[j], col, pivot_idx, k);
                }
                delete_column(col);
            }
        }
    }
    double end_time = MPI_Wtime();
    if (gj.my_rank == 0)
        printf("%f\n", end_time - start_time);
    count_processes();
}

/**
 * Checks if the column given belongs to the MPI process passed
 * @param col_idx The column idx
 * @return On true, 1 is returned, otherwise 0
 */
int is_my_column(int col_idx) {
    return process_of_column(col_idx) == gj.my_rank;
}

/**
 * Returns the ID (rank) of the MPI Process that owns the given column id
 * @param col_idx The column id (from 0..n)
 * @return The rank of the MPI Process
 */
int process_of_column(int col_idx) {
    return (gj.use_group_distribution) ?
        col_idx / gj.group_number : col_idx % gj.proc_num;
}

/**
 * Returns the number of columns that a process has
 * (group number := dimension/process_number)
 * Rank_0 has group_number + 1, all others group_number
 * @return The number of columns of the rank-th process
 */
int number_of_cols() {
    return (gj.my_rank == 0) ? gj.group_number + 1 : gj.group_number;
}

/**
 * Maps the global idx (k - step) to the idx of the local array
 * @param k The k-step idx
 * @return The local idx
 */
int map_global_to_local(int k) {
    return ((gj.use_group_distribution) ? k % gj.group_number : k / gj.proc_num);
}

/**
 * Returns the number of recipients depending on the distribution
 * @return The number of recipients
 */
int number_of_recipients() {
    int recipients_num = (gj.use_group_distribution) ? gj.proc_num - gj.my_rank : gj.proc_num - 1;
    // If it's the first process (non master), it must send to [proc_num - rank]
    // processes
    if (gj.use_group_distribution && gj.my_rank == 0)
        recipients_num++;
    return recipients_num;
}

/**
 * Send the k-th column to all the > k column owners
 * @param k The current step of the method
 * @param max_idx The pivot element index
 * @param my_cols The columns array of the current MPI process
 * @return Array of MPI_Requests, that was populated on Isend
 */
MPI_Request* send_column(int k, int max_idx, column_t* my_cols) {
    // Get the k-th column
    int local_idx = map_global_to_local(k);
    // Serialize the struct to the local array
    gj.dummy_col[0] = max_idx;
    memcpy(&(gj.dummy_col[1]), my_cols[local_idx]->data, gj.dimension * sizeof(float));
    // Send it to all the > k owners, but gather indices to my own
    // (that are > k)
    MPI_Request *requests = malloc(
            number_of_recipients() * sizeof(MPI_Request));
    int rank, requests_idx = 0;
    if (gj.use_group_distribution) {
        for (rank = gj.my_rank + 1; rank < gj.proc_num; rank++) {
            MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, rank, k,
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
        }
        if (gj.my_rank != 0) {
            // Send to master process as well
            MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, 0, k,
                        MPI_COMM_WORLD, &(requests[requests_idx++]));
        }
    }
    else {

        if (k >= gj.dimension - gj.proc_num + 2) {
            for (rank = gj.my_rank + 1; rank < gj.proc_num; rank++)
                MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, rank, k,
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
            // Send to master process
            MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, 0, k,
                    MPI_COMM_WORLD, &(requests[requests_idx++]));
        }
        else {
            for (rank = 0; rank < gj.proc_num; rank++) {
                if (rank != gj.my_rank)
                    MPI_Isend(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, rank, k,
                        MPI_COMM_WORLD, &(requests[requests_idx++]));
            }
        }
    }
    gj.recipients_number = requests_idx;
    return requests;
}

/**
 * Receives the column k and saves it to the dummy_col field
 * @param k The current step of the Gauss Jordan method
 * @return The newly arrived column
 */
column_t receive_column(int k) {
    MPI_Request request;
    MPI_Irecv(gj.dummy_col, gj.dimension + 1, MPI_FLOAT, process_of_column(k),
            k, MPI_COMM_WORLD, &request);
    wait_wrapper(&request);
    return create_column(k, &(gj.dummy_col[1]));
}

/**
 * Randomly creates a [n x (n + 1)] matrix
 * @param dimension The dimension n of the matrix
 * @return The augmented matrix
 */
float** create_augmented_matrix(int dimension) {
    int i,j;
    float** matrix = malloc(dimension * sizeof(float*));
    for (i = 0; i < dimension; i++) {
        matrix[i] = malloc((dimension + 1) * sizeof(float));
        for (j = 0; j < dimension + 1; j++) {
            matrix[i][j] = rand();
        }
    }
    return matrix;
}

void destroy_augmented_matrix(float** augmented_m) {
    int i;
    for (i = 0; i < gj.dimension; i++) {
        free(augmented_m[i]);
        augmented_m[i] = NULL;
    }
    free(augmented_m);
    augmented_m = NULL;
}

void count_processes() {
    int local_int = 1, total_int = 0;
    MPI_Reduce(&local_int, &total_int, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (gj.my_rank == 0) {
        printf("Master process reporting [%d] number of processes\n", total_int);
        if (total_int != gj.proc_num) {
            printf("ERROR:\tproc_num: %d\n", gj.proc_num);
        }
    }
}
