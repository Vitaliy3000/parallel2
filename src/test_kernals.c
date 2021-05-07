#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "kernals.h"


const double EPS = 0.00001;
const int COUNT_CALLS = 10;
const int SIZE_OF_ARRAY = 10000000;
const int MAX_COUNT_THREADS = 13;


int compare(double x, double y) {
    return ( (x > y - EPS) && (x < y + EPS) );
}


int test_sum() {
    double x[] = {1, 2, 3};
    return compare(sum(3, x), 6);
}


double mean_time_sum() {
    double *x = malloc(SIZE_OF_ARRAY * sizeof(double));
    for (int i = 0; i < SIZE_OF_ARRAY; i++) x[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) sum(SIZE_OF_ARRAY, x);
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int test_L2() {
    double x[] = {1, 2, 3};
    return compare(L2(3, x), 3.74165738677);
}


double mean_time_L2() {
    double *x = malloc(SIZE_OF_ARRAY * sizeof(double));
    for (int i = 0; i < SIZE_OF_ARRAY; i++) x[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) L2(SIZE_OF_ARRAY, x);
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int test_dot() {
    double x[] = {1, 2, 3};
    double y[] = {2, 3, 4};
    return compare(dot(3, x, y), 20);
}


double mean_time_dot() {
    double *x = malloc(SIZE_OF_ARRAY * sizeof(double));
    double *y = malloc(SIZE_OF_ARRAY * sizeof(double));
    for (int i = 0; i < SIZE_OF_ARRAY; i++) x[i] = i / 10000.0;
    for (int i = 0; i < SIZE_OF_ARRAY; i++) y[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) dot(SIZE_OF_ARRAY, x, y);
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int test_axpby() {
    double x[] = {1, 2, 3};
    double y[] = {2, 3, 4};
    double res[] = {-7, -9, -11};
    axpby(3, 3, x, -5, y);
    return ( (compare(sum(3, x), sum(3, res))) && (compare(L2(3, x), L2(3, res))) );
}


double mean_time_axpby() {
    double *x = malloc(SIZE_OF_ARRAY * sizeof(double));
    double *y = malloc(SIZE_OF_ARRAY * sizeof(double));
    for (int i = 0; i < SIZE_OF_ARRAY; i++) x[i] = i / 10000.0;
    for (int i = 0; i < SIZE_OF_ARRAY; i++) y[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) axpby(SIZE_OF_ARRAY, 3, x, -5, y);
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int test_VVbe() {
    double x[] = {1, 2, 3};
    double y[] = {2, 3, 4};
    double res[] = {2, 6, 12};
    double z[3];
    VVbe(3, x, y, z);
    return ( (compare(sum(3, z), sum(3, res))) && (compare(L2(3, z), L2(3, res))) );
}


double mean_time_VVbe() {
    double *x = malloc(SIZE_OF_ARRAY * sizeof(double));
    double *y= malloc(SIZE_OF_ARRAY * sizeof(double));
    double *z = malloc(SIZE_OF_ARRAY * sizeof(double));
    for (int i = 0; i < SIZE_OF_ARRAY; i++) x[i] = i / 10000.0;
    for (int i = 0; i < SIZE_OF_ARRAY; i++) y[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) VVbe(SIZE_OF_ARRAY, x, y, z);
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int test_SpMV() {
    int IA[] = {0, 2, 4, 5, 6, 10, 13};
    int JA[] = {0, 3, 1, 5, 0, 3, 0, 1, 2, 3, 0, 4, 5};
    double A[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    double b[] = {1, 2, 4, 8, 1, 3};
    double c[13];
    double res[] = {17, 18, 5, 48, 139, 62};
    SpMV(6, IA, JA, A, b, c);
    return ( (compare(sum(3, c), sum(3, res))) && (compare(L2(3, c), L2(3, res))) );
}


double mean_time_SpMV() {
    int *IA = malloc((SIZE_OF_ARRAY + 1) * sizeof(int));
    int *JA = malloc(SIZE_OF_ARRAY * sizeof(int));
    double *A = malloc(SIZE_OF_ARRAY * sizeof(double));;
    double *b = malloc(SIZE_OF_ARRAY * sizeof(double));;
    double *c = malloc(SIZE_OF_ARRAY * sizeof(double));;
    for (int i = 0; i < SIZE_OF_ARRAY+1; i++) IA[i] = i;
    for (int i = 0; i < SIZE_OF_ARRAY; i++) JA[i] = i;
    for (int i = 0; i < SIZE_OF_ARRAY; i++) A[i] = i / 10000.0;
    for (int i = 0; i < SIZE_OF_ARRAY; i++) b[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) SpMV(SIZE_OF_ARRAY, IA, JA, A, b, c);;
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int test_copy() {
    double x[] = {1, 2, 3};
    double y[3];
    double res[] = {1, 2, 3};
    copy(3, x, y);
    return ( (compare(sum(3, y), sum(3, res))) && (compare(L2(3, y), L2(3, res))) );
}


double mean_time_copy() {
    double *x = malloc(SIZE_OF_ARRAY * sizeof(double));
    double *y = malloc(SIZE_OF_ARRAY * sizeof(double));
    for (int i = 0; i < SIZE_OF_ARRAY; i++) x[i] = i / 10000.0;
    double start_time = omp_get_wtime();
    for (int i = 0; i < COUNT_CALLS; i++) copy(SIZE_OF_ARRAY, x, y);
    return (omp_get_wtime() - start_time) / COUNT_CALLS;
}


int main(int argc, char *argv[]) {
    printf("\n");

    assert(test_sum());
    printf("Test sum: Mean time: %.4f\n", mean_time_sum());

    assert(test_L2());
    printf("Test L2: Mean time: %.4f\n", mean_time_L2());

    assert(test_dot());
    printf("Test dot: Mean time: %.4f\n", mean_time_dot());

    assert(test_axpby());
    printf("Test axpby: Mean time: %.4f\n", mean_time_axpby());

    assert(test_VVbe());
    printf("Test VVbe: Mean time: %.4f\n", mean_time_VVbe());

    assert(test_SpMV());
    printf("Test SpMV: Mean time: %.4f\n", mean_time_SpMV());

    assert(test_copy());
    printf("Test copy: Mean time: %.4f\n", mean_time_copy());

    printf("\n");

    return 0;
}
