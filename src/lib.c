#include "lib.h"
#include <math.h>
#include <omp.h>


double sum(int N, double x[]) {
    double result = 0;
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < N; i++) result += x[i];
    return result;
}


double L2(int N, double x[]) {
    return sqrtf(dot(N, x, x));
}


double dot(int N, double x[], double y[]) {
    double result = 0;
    #pragma omp parallel for reduction(+:result)
    for (int i = 0; i < N; i++) result += x[i] * y[i];
    return result;
}


void axpby(int N, double a, double x[], double b, double y[]) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) x[i] = a * x[i] + b * y[i];
}


void VVbe(int N, double x[], double y[], double z[]) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) z[i] = x[i] * y[i];
}


void SpMV(int N, int IA[], int JA[], double A[], double b[], double c[]) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        c[i] = 0;
        for (int j=IA[i]; j < IA[i+1]; j++)
            c[i] += A[j] * b[JA[j]];
    }
}


void copy(int N, double x[], double y[]) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) y[i] = x[i];
}
