
#ifndef COLUMN_H
#define	COLUMN_H

#ifdef	__cplusplus
extern "C" {
#endif

    struct Column_Type {
        float* data;
        int idx;
    };
    typedef struct Column_Type* column_t;

    column_t create_column(int idx, float* data);
    void delete_column(column_t col_ptr);
    void print_column(column_t column);
    int pivot(column_t column, int k);
    void modify(column_t col, column_t pivot_col,
    int pivot_idx, int k);


#ifdef	__cplusplus
}
#endif

#endif	/* COLUMN_H */

