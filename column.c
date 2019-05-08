#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "column.h"
#include "gauss_jordan.h"

extern gauss_jordan gj;

// Forward declarations
void swap_elements(column_t column, int i, int j);

/**
 * Creates a Column datatype (allocates the memory
 * needed and returns a pointer to it)
 * @param idx The index of the column in the original matrix
 * @param data The coefficients of the column
 * @return A pointer to the allocated Column datatype
 */
column_t create_column(int idx, float* data) {
    column_t col = malloc(sizeof(struct Column_Type));
    col->idx = idx;
    col->data = malloc(gj.dimension * sizeof(float));
    memcpy(col->data, data, gj.dimension * sizeof(float));
    return col;
}

/**
 * Deletes the space allocated for the column and nullifies the
 * pointer to it
 * @param col_ptr The address of the pointer to the column datatype
 */
void delete_column(column_t col_ptr) {
    (col_ptr)->idx = -1;
    free((col_ptr)->data);
    (col_ptr)->data = NULL;
}

/**
 * Prints the column of length dimension (defined by the global 
 * Gauss Jordan struct)
 * @param column A pointer to the column to be printed
 */
void print_column(column_t column) {
    if (column->data != NULL) {
        printf("Column[%d]: ", column->idx);
        int i;
        for (i = 0; i < gj.dimension; i++) {
            printf("%f | ", column->data[i]);
        }
        printf("\n");
    }
    else {
        printf("Column isn't instantiated\n");
    }
}

/**
 * Finds max element from k to n, swaps
 * the max element with the [k,k] element
 * @param column The column to be processed
 * @param k The current step of the Gauss Jordan method
 * @return The idx of max (pivot) element 
 */
int pivot(column_t column, int k) {
    // Find max element
    float max = column->data[k];
    int max_idx = k;
    int i;
    for (i = k + 1; i < gj.dimension; i++) {
        if (column->data[i] > max) {
            max = column->data[i];
            max_idx = i;
        }
    }
    // Swap the max element with the [k,k] element
    swap_elements(column, max_idx, k);
    return max_idx;
}

/**
 * Plain old swap method
 * @param column The column which elements are to be swapped
 * @param i The idx of the 1st element
 * @param j The idx of the 2nd element
 */
void swap_elements(column_t column, int i, int j) {
    float tmp;
    tmp = column->data[i];
    column->data[i] = column->data[j];
    column->data[j] = tmp;
}

/**
 * Modifies the contents of the given column using the
 * k-column (the pivot one)
 * @param column The column to be modified
 * @param pivot_col The k-column of the matrix
 * @param pivot_idx The idx of the pivot element in the k-column
 */
void modify(column_t col, column_t pivot_col, 
        int pivot_idx, int k) {
    swap_elements(col, k, pivot_idx);
    col->data[k] = col->data[k]/pivot_col->data[k];
    int i;
    for (i = 0; i < gj.dimension; i++) {
        if (i != k)
            col->data[i] = 
                    col->data[i] - pivot_col->data[i] * col->data[k];
    }
}