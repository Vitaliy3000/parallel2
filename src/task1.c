#include <stdlib.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "parser.h"
#include "lib.h"
#include "matrix.h"
#include "kernals_omp.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))


int DEBUG = 0;

const int LOW_K1 = 0;
const int HIGH_K1 = 10000;
const int LOW_K2 = 0;
const int HIGH_K2 = 10000;
const int LOW_Nx = 1;
const int HIGH_Nx = 10000;
const int LOW_Ny = 1;
const int HIGH_Ny = 10000;
const int LOW_maxiter = 1;
const int HIGH_maxiter = 10000;
const double LOW_tol = 0;
const double HIGH_tol = 1;
const int LOW_count_threads = 1;
const int HIGH_count_threads = 100;


char* HELP_STRING = "\
  usage:\n\
    task1 [--help]\n\
    task1 --Nx <Nx> --Ny <Ny> --K1 <K1> --K2 <K2> -n <n> --maxiter <N> --tol <tol> [--file  <filename>]\n\
  options:\n\
    [--help]     Show this screen.\n\
    --Nx         Count of elements in horizontal\n\
    --Ny         Count of elements in vertical\n\
    --K1         Count of squares\n\
    --K2         Count of cut squares\n\
    -n           Count of threads\n\
    --maxiter    Maximum of count steps in solver\n\
    --tol        Relative tolerance for termination.\n\
    [--filename] Debug mode\n\
";


struct res_by_kernel {
    double time_dot;
    double time_axpby;
    double time_SpMV;
    double time_VVbe;

    double gflops_dot;
    double gflops_axpby;
    double gflops_SpMV;
    double gflops_VVbe;
};


void parse_args(int argc, char *argv[], int *Nx, int *Ny, int *K1, int *K2, int *maxiter, double *tol, int *count_threads, char *filename) {
    if ((argc == 1) || (strcmp(argv[1], "--help"))==0) {
        printf(HELP_STRING);
        exit(0);
    }

    parse_int_arg(argc, argv, "--Nx", Nx, LOW_Nx, HIGH_Nx);
    parse_int_arg(argc, argv, "--Ny", Ny, LOW_Ny, HIGH_Ny);
    parse_int_arg(argc, argv, "--K1", K1, LOW_K1, HIGH_K1);
    parse_int_arg(argc, argv, "--K2", K2, LOW_K2, HIGH_K2);
    parse_int_arg(argc, argv, "--maxiter", maxiter, LOW_maxiter, HIGH_maxiter);
    parse_int_arg(argc, argv, "-n", count_threads, LOW_count_threads, HIGH_count_threads);
    parse_double_arg(argc, argv, "--tol", tol, LOW_tol, HIGH_tol);

    parse_str_arg(argc, argv, "--filename", filename);
}


void calc_count_figures(int Nx, int Ny, int K1, int K2, int *count_triangles, int *count_squares) {
    *count_triangles = 2 * ( K2 * ( Nx * Ny / (K1 + K2) ) + MAX(0, Nx * Ny % (K1 + K2) - K1) );
    *count_squares = K1 * ( Nx * Ny / (K1 + K2) ) + MIN(K1, Nx * Ny % (K1 + K2));
}


void gen_graph(int K, int *offset_elements, int *elements, int Nx, int Ny, int K1, int K2) {
    int SWITCH_FLAG = 0, current_x = 0, current_y = 0, k = 0, p = 0;

    offset_elements[0] = 0;

    while (k < K) {
        if (SWITCH_FLAG >= K1) {
            elements[p++] = current_y*(Nx+1) + current_x;
            elements[p++] = current_y*(Nx+1) + current_x + 1;
            elements[p++] = (current_y+1)*(Nx+1) + current_x;

            offset_elements[k+1] = offset_elements[k] + 3;
            ++k;

            elements[p++] = current_y*(Nx+1) + current_x + 1;
            elements[p++] = (current_y+1)*(Nx+1) + current_x + 1;
            elements[p++] = (current_y+1)*(Nx+1) + current_x;

            offset_elements[k+1] = offset_elements[k] + 3;
            ++k;
        } else {
            elements[p++] = current_y*(Nx+1) + current_x;
            elements[p++] = current_y*(Nx+1) + current_x + 1;
            elements[p++] = (current_y+1)*(Nx+1) + current_x + 1;
            elements[p++] = (current_y+1)*(Nx+1) + current_x;
            offset_elements[k+1] = offset_elements[k] + 4;
            ++k;
        }

        ++current_x;

        if (current_x == Nx) {
            ++current_y;
            current_x = 0;
        }

        ++SWITCH_FLAG;

        if (SWITCH_FLAG == (K1 + K2)) SWITCH_FLAG = 0;
    }
}


void solve(int N, int IA[], double x[], int JA[], double A[], double b[], int *k, double *res, int maxiter, double tol, struct res_by_kernel *gflops) {
    /*
    N - length of IA, x, b
    n - count of iter
    res - norm L2
    */
    double start_time;

    gflops->gflops_dot = 0;
    gflops->gflops_axpby = 0;
    gflops->gflops_SpMV = 0;
    gflops->gflops_VVbe = 0;

    gflops->time_dot = 0;
    gflops->time_axpby = 0;
    gflops->time_SpMV = 0;
    gflops->time_VVbe = 0;

    if (DEBUG) printf("\nDEBUG solve:\n");

    double alpha, betta, rho, last_rho = 1;
    double *r = (double*) strict_calloc(N, sizeof(double));
    double *z = (double*) strict_calloc(N, sizeof(double));
    double *p = (double*) strict_malloc(N * sizeof(double));
    double *q = (double*) strict_malloc(N * sizeof(double));
    double *M = (double*) strict_malloc(N * sizeof(double));

    for (int i = 0; i < N; i++)
        for (int j = IA[i]; j < IA[i+1]; j++)
            if (i == JA[j])
                M[i] = 1/A[j];

    start_time = omp_get_wtime();
    SpMV(N, IA, JA, A, x, r);
    gflops->time_SpMV += omp_get_wtime() - start_time;
    gflops->gflops_SpMV += 2 * 1E-9 * IA[N];

    start_time = omp_get_wtime();
    axpby(N, -1, r, 1, b);
    gflops->time_axpby += omp_get_wtime() - start_time;
    gflops->gflops_axpby += 3 * 1E-9 * N;

    *k = 1;
    while (1) {
        start_time = omp_get_wtime();
        VVbe(N, M, r, z);
        gflops->time_VVbe += omp_get_wtime() - start_time;
        gflops->gflops_VVbe += 1E-9 * N;

        start_time = omp_get_wtime();
        rho = dot(N, r, z);
        gflops->time_dot += omp_get_wtime() - start_time;
        gflops->gflops_dot += 2 * 1E-9 * N;

        if (*k == 1) {
            copy(N, z, p);
        } else {
            betta = rho / last_rho;

            start_time = omp_get_wtime();
            axpby(N, betta, p, 1, z);
            gflops->time_axpby += omp_get_wtime() - start_time;
            gflops->gflops_axpby += 3 * 1E-9 * N;
        }

        start_time = omp_get_wtime();
        SpMV(N, IA, JA, A, p, q);
        gflops->time_SpMV += omp_get_wtime() - start_time;
        gflops->gflops_SpMV += 2 * 1E-9 * IA[N];

        start_time = omp_get_wtime();
        alpha = rho / dot(N, p, q);
        gflops->time_dot += omp_get_wtime() - start_time;
        gflops->gflops_dot += 2 * 1E-9 * N;

        start_time = omp_get_wtime();
        axpby(N, 1, x, alpha, p);
        gflops->time_axpby += omp_get_wtime() - start_time;
        gflops->gflops_axpby += 3 * 1E-9 * N;

        start_time = omp_get_wtime();
        axpby(N, 1, r, -alpha, q);
        gflops->time_axpby += omp_get_wtime() - start_time;
        gflops->gflops_axpby += 3 * 1E-9 * N;

        *res = L2(N, r);

        printf("\nNumber of iteration: %d\tres = %.5f\n", *k, *res);

        if ((*res / L2(N, b) < tol) || (*k >= maxiter))
            break;
        else
            ++(*k);

        last_rho = rho;
    }

    free(r);
    free(z);
    free(p);
    free(q);
    free(M);
}


int main(int argc, char *argv[]) {
    int k, Nx, Ny, K1, K2, K, count_triangles, count_squares, exepted_memory=0, K_full_size, flag = 0, count_steps, maxiter, count_threads;
    double sum, tol;
    int *elements, *offset_elements;
    double res;
    FILE *f;
    char filename[MAX_SIZE_OF_STRING] = "";
    double start_time, end_time;
    struct res_by_kernel gflops;

    start_time = omp_get_wtime();

    parse_args(argc, argv, &Nx, &Ny, &K1, &K2, &maxiter, &tol, &count_threads, filename);

    if (ERROR_FLAG) exit(EXIT_FAILURE);

    if (filename[0] != '\0') {
        f = fopen(filename, "w");
        DEBUG = 1;
    }

    if (DEBUG) {
        fprintf(f, "Nx = %d\n", Nx);
        fprintf(f, "Ny = %d\n", Ny);
        fprintf(f, "K1 = %d\n", K1);
        fprintf(f, "K2 = %d\n", K2);
        fprintf(f, "maxiter = %d\n", maxiter);
        fprintf(f, "tol = %.5f\n", tol);
        fprintf(f, "filename = %s\n", filename);
    }

    calc_count_figures(Nx, Ny, K1, K2, &count_triangles, &count_squares);
    K = count_squares + count_triangles;
    K_full_size = (4*count_squares + 3*count_triangles);

    offset_elements = (int*) strict_malloc(sizeof(int) * (K+1));
    exepted_memory += sizeof(int) * (K+1);

    elements = (int*) strict_malloc(sizeof(int) * K_full_size);
    exepted_memory += sizeof(int) * K_full_size;

    gen_graph(K, offset_elements, elements, Nx, Ny, K1, K2);

    int count_edges = Ny*(Nx+1) + Nx*(Ny+1) + count_triangles / 2;
    int count_nodes = (Nx+1) * (Ny+1);

    if (DEBUG) {
        fprintf(f, "\nGraph:\n");
        for (int k = 0; k < K; k++) {
            fprintf(f, "    Element:");
            for (int i = offset_elements[k]; i < offset_elements[k+1]; ++i)
                fprintf(f, " %d", elements[i]);

            fprintf(f, "\n");
        }

        fprintf(f, "\n");

        fprintf(f, "    Count nodes: %d\n", count_nodes);
        fprintf(f, "    Count edges: %d\n", count_edges);

        int *array = (int*) strict_calloc(count_nodes * count_nodes, sizeof(int));

        graph2matrix(K, offset_elements, elements, count_nodes, array);

        fprintf(f, "\nGraph2Matrix:\n");
        for (int i = 0; i < count_nodes; i++) {
            for (int j = 0; j < count_nodes; j++)
                fprintf(f, "%d ", array[i*count_nodes + j]);
            fprintf(f, "\n");
        }
    }

    int *IA = (int*) strict_malloc((count_nodes + 1) * sizeof(int));
    exepted_memory += (count_nodes + 1) * sizeof(int);
    int *JA = (int*) strict_malloc((2 * count_edges + count_nodes) * sizeof(int));
    exepted_memory += (2 * count_edges + count_nodes) * sizeof(int);

    graph2csr(K, offset_elements, elements, count_edges, count_nodes, IA, JA);

    if (DEBUG) {
        fprintf(f, "\nCSR:");
        fprintf(f, "\n    IA:");
        for (int i = 0; i < count_nodes + 1; i++)
            fprintf(f, " %d", IA[i]);
        fprintf(f, "\n    JA:");
        for (int i = 0; i < 2 * count_edges + count_nodes; i++)
            fprintf(f, " %d", JA[i]);

        int *array = (int*) strict_calloc(count_nodes * count_nodes, sizeof(int));

        csr2matrix(count_nodes + 1, IA, 2 * count_edges + count_nodes, JA, array);

        fprintf(f, "\n\nCSR2Matrix:\n");
        for (int i = 0; i < count_nodes; i++) {
            for (int j = 0; j < count_nodes; j++)
                fprintf(f, "%d ", array[i*count_nodes + j]);
            fprintf(f, "\n");
        }
    }

    end_time = omp_get_wtime();

    printf("\nGraph in CSR format was build\n");
    printf("Time of building: %.5f seconds\n", end_time - start_time);
    printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

    free(elements);
    exepted_memory -= sizeof(int) * K_full_size;
    free(offset_elements);
    exepted_memory -= sizeof(int) * (K+1);

    start_time = omp_get_wtime();

    double *A = (double*) strict_malloc((2 * count_edges + count_nodes) * sizeof(double));
    exepted_memory += (2 * count_edges + count_nodes) * sizeof(double);

    k = 0;
    for (int i = 0; i < count_nodes; i++) {
        sum = 0;
        for (int j = IA[i]; j < IA[i+1]; j++) {
            if (i != JA[j]) {
                A[k] = cos(i*JA[j] + i + JA[j]);
                sum += fabs(A[k]);
            } else {
                flag = k;
            }
            ++k;
        }
        A[flag] = 1.234*sum;
    }

    double *b = (double*) strict_malloc(count_nodes * sizeof(double));
    exepted_memory += count_nodes * sizeof(double);

    for (int i = 0; i < count_nodes; i++) b[i] = sin(i);

    end_time = omp_get_wtime();

    printf("\nCalculated values in matrix\n");
    printf("Time of calculating: %.5f seconds\n", end_time - start_time);
    printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

    start_time = omp_get_wtime();

    double *x = (double*) strict_calloc(count_nodes, sizeof(double));
    exepted_memory += count_nodes * sizeof(double);

    omp_set_num_threads(count_threads);

    solve(count_nodes, IA, x, JA, A, b, &count_steps, &res, maxiter, tol, &gflops);

    if (DEBUG) {
        fprintf(f, "\n\nResults:\n");
        fprintf(f, "  res: %.5f\n", res);
        fprintf(f, "  count steps: %d\n", count_steps);
        fprintf(f, "  X: ");
        for (int i = 0; i < count_nodes; i++) fprintf(f, " %.4f", x[i]);
    }

    end_time = omp_get_wtime();

    printf("\nSolved\n");
    printf("Time of calculating: %.5f seconds\n", end_time - start_time);
    printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

    printf("\nTiming by kernels\n");
    printf("Time of dot: %.5f seconds, gflops: %.5f\n", gflops.time_dot, gflops.gflops_dot/gflops.time_dot);
    printf("Time of axpby: %.5f seconds, gflops: %.5f\n", gflops.time_axpby, gflops.gflops_axpby/gflops.time_axpby);
    printf("Time of SpMV: %.5f seconds, gflops: %.5f\n", gflops.time_SpMV, gflops.gflops_SpMV/gflops.time_SpMV);
    printf("Time of VVbe: %.5f seconds, gflops: %.5f\n\n", gflops.time_VVbe, gflops.gflops_VVbe/gflops.time_VVbe);
    printf(
        "Time of solver: %.5f seconds, gflops: %.5f\n\n",
        gflops.time_dot+gflops.time_axpby+gflops.time_SpMV+gflops.time_VVbe,
        (gflops.gflops_dot+gflops.gflops_axpby+gflops.gflops_SpMV+gflops.gflops_VVbe)
        / (gflops.time_dot+gflops.time_axpby+gflops.time_SpMV+gflops.time_VVbe)
    );

    return 0;
}
