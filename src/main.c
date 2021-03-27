#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "lib.h"

#define MAX_SIZE_OF_STRING 40
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))


int DEBUG = 0;

char* MAX_SIZE_OF_STRING_ERROR = "Max size of argumet %s must be less 40 letters\n";
char* DUPLICATE_ARGUMENT_ERROR = "Duplicate argument of %s\n";
char* NOT_FOUND_ARGUMENT_ERROR = "Not found argument of %s\n";
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

char* LOW_ARGUMENT_ERROR_INTEGER = "Argument %s must be bigger %d\n";
char* HIGH_ARGUMENT_ERROR_INTEGER = "Argument %s must be lower %d\n";
char* LOW_ARGUMENT_ERROR_FLOAT = "Argument %s must be bigger %.2f\n";
char* HIGH_ARGUMENT_ERROR_FLOAT = "Argument %s must be lower %.2f\n";
char* INTEGER_ARGUMENT_ERROR = "Argument %s must be ineger\n";
char* FLOAT_ARGUMENT_ERROR = "Argument %s must be double\n";

char* HELP_STRING = "\
  usage:\n\
    a.exe [--help]\n\
    a.exe --Nx <Nx> --Ny <Ny> --K1 <K1> --K2 <K2> -n <n> --maxiter <N> --tol <tol> [--file  <filename>]\n\
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

int ERROR_FLAG = 0;


struct time_by_kernel {
    double dot;
    double axpby;
    double SpMV;
    double VVbe;
};


void* strict_malloc(size_t size) {
    void *new_ptr = malloc( size);
 
    if (new_ptr == NULL) {
        perror("calloc return NULL");
        exit(EXIT_FAILURE);
    }

    return new_ptr;
}


void* strict_calloc(size_t num, size_t size) {
    void *new_ptr = calloc(num, size);
 
    if (new_ptr == NULL) {
        perror("calloc return NULL");
        exit(EXIT_FAILURE);
    }

    return new_ptr;
}


void* strict_realloc(void *ptr, size_t newsize) {
    void *new_ptr = realloc(ptr, newsize);

    if (new_ptr == NULL) {
        perror("realloc return NULL");
        exit(EXIT_FAILURE);
    }

    return new_ptr;
}


int strict_atoi(char str[], int *ptr) {
    for(int i = 0; str[i] != '\0'; i++) {
        if (((str[i] > '9') || (str[i] < '0')) && (str[i] != '-')) {
            return -1;
        }
    }

    *ptr = atoi(str);
    return 0;
}


double strict_atof(char str[], double *ptr) {
    for(int i = 0; str[i] != '\0'; i++) {
        if (((str[i] > '9') || (str[i] < '0')) && (str[i] != '-') && (str[i] != '.')) {
            return -1;
        }
    }

    *ptr = atof(str);
    return 0;
}


void parse_int_arg(int argc, char *argv[], char expected_key[], int ptr[], int low, int high) {
    char value[MAX_SIZE_OF_STRING] = "";
    
    for (int i = 1; i < argc-1; i++) {
        if (strcmp(argv[i], expected_key) == 0) {
            if (value[0] != '\0') {
                fprintf(stderr, DUPLICATE_ARGUMENT_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            if (strlen(argv[i+1]) >= MAX_SIZE_OF_STRING) {
                fprintf(stderr, MAX_SIZE_OF_STRING_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            strcpy(value, argv[i+1]);
        }
    }

    if (value[0] == '\0') {
        fprintf(stderr, NOT_FOUND_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (strict_atoi(value, ptr) == -1) {
        fprintf(stderr, INTEGER_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr < low) {
        fprintf(stderr, LOW_ARGUMENT_ERROR_INTEGER, expected_key, low);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr > high) {
        fprintf(stderr, HIGH_ARGUMENT_ERROR_INTEGER, expected_key, high);
        ERROR_FLAG = 1;
        return;
    }

    return;
}


void parse_double_arg(int argc, char *argv[], char expected_key[], double ptr[], double low, double high) {
    char value[MAX_SIZE_OF_STRING] = "";

    for (int i = 1; i < argc-1; i++) {
        if (strcmp(argv[i], expected_key) == 0) {
            if (value[0] != '\0') {
                fprintf(stderr, DUPLICATE_ARGUMENT_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            if (strlen(argv[i+1]) >= MAX_SIZE_OF_STRING) {
                fprintf(stderr, MAX_SIZE_OF_STRING_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            strcpy(value, argv[i+1]);
        }
    }

    if (value[0] == '\0') {
        fprintf(stderr, NOT_FOUND_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (strict_atof(value, ptr) == -1) {
        fprintf(stderr, FLOAT_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr < low) {
        fprintf(stderr, LOW_ARGUMENT_ERROR_FLOAT, expected_key, low);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr > high) {
        fprintf(stderr, HIGH_ARGUMENT_ERROR_FLOAT, expected_key, high);
        ERROR_FLAG = 1;
        return;
    }

    return;
}


void parse_str_arg(int argc, char *argv[], char expected_key[], char *str) {
    for (int i = 1; i < argc-1; i++) {
        if (strcmp(argv[i], expected_key) == 0) {
            if (str[0] != '\0') {
                fprintf(stderr, DUPLICATE_ARGUMENT_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            if (strlen(argv[i+1]) >= MAX_SIZE_OF_STRING) {
                fprintf(stderr, MAX_SIZE_OF_STRING_ERROR, expected_key);
                ERROR_FLAG = 1;
                return;
            }

            strcpy(str, argv[i+1]);
        }
    }
}


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


int cmpfunc(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}


void graph2csr(int K, int offset_elements[], int elements[], int count_edges, int count_nodes, int IA[],  int JA[]) {
    int left, right, flag;

    int *counter_IA_d = strict_calloc(count_nodes, sizeof(int));

    for (int k = 0; k < K; k++) {
        int i = offset_elements[k];
        for (; i < offset_elements[k+1] - 1; ++i) {
            left = elements[i];
            right = elements[i+1];

            ++counter_IA_d[left];
            ++counter_IA_d[right];
        }

        left = elements[i];
        right = elements[offset_elements[k]];

        ++counter_IA_d[left];
        ++counter_IA_d[right];
    }

    int **NE_d = strict_malloc(count_nodes * sizeof(int*));
    for (int i = 0; i < count_nodes; i++) NE_d[i] = strict_malloc(counter_IA_d[i] * sizeof(int));
    for (int i = 0; i < count_nodes; i++) counter_IA_d[i] = 0;

    for (int k = 0; k < K; k++) {
        int i = offset_elements[k];
        for (; i < offset_elements[k+1] - 1; ++i) {
            left = elements[i];
            right = elements[i+1];

            NE_d[left][counter_IA_d[left]] = right;
            NE_d[right][counter_IA_d[right]] = left;

            ++counter_IA_d[left];
            ++counter_IA_d[right];
        }

        left = elements[i];
        right = elements[offset_elements[k]];

        NE_d[left][counter_IA_d[left]] = right;
        NE_d[right][counter_IA_d[right]] = left;

        ++counter_IA_d[left];
        ++counter_IA_d[right];
    }


    for (int i = 0; i < count_nodes; i++) qsort(NE_d[i], counter_IA_d[i], sizeof(int), cmpfunc);

    int *counter_IA = strict_calloc(count_nodes, sizeof(int));

    for (int i = 0; i < count_nodes; i++) {
        for (int j = 0; j < counter_IA_d[i]-1; j++) {
            if (NE_d[i][j] != NE_d[i][j+1]) {
                ++counter_IA[i];
            }
        }

        ++counter_IA[i];
    }

    int **NE = strict_malloc(count_nodes * sizeof(int*));
    for (int i = 0; i < count_nodes; i++) NE[i] = strict_malloc(counter_IA[i] * sizeof(int));
    for (int i = 0; i < count_nodes; i++) counter_IA[i] = 0;

    for (int i = 0; i < count_nodes; i++) {
        int j = 0;
        for (; j < counter_IA_d[i]-1; j++) {
            if (NE_d[i][j] != NE_d[i][j+1]) {
                NE[i][counter_IA[i]] = NE_d[i][j];
                ++counter_IA[i];
            }
        }

        NE[i][counter_IA[i]] = NE_d[i][j];
        ++counter_IA[i];
    }

    IA[0] = 0;
    for (int i = 1; i < count_nodes+1; i++) IA[i] = IA[i-1] + counter_IA[i-1] + 1;
    for (int i = 0; i < count_nodes; i++) counter_IA[i] = 0;


    int k = 0;
    for (int i = 0; i < count_nodes; i++) {
        flag = 0;
        for (int j = 0; j < IA[i+1] - IA[i]; j++) {
            if ( 
                (flag == 0)
                && (
                    (j + 1 == IA[i+1] - IA[i])
                    || ( (j == 0) && (i < NE[i][j]) )
                    || ( (i > NE[i][j]) || (i < NE[i][j+1]) ) 
                ) 
            ) {
                JA[k] = i;
                ++k;
                ++j;
                flag = 1;
            }
            JA[k] = NE[i][j-flag];
            ++k;
        }
    }

    free(counter_IA_d);
    free(counter_IA);

    for (int i = 0; i < count_nodes; i++) free(NE_d[i]);
    free(NE_d);

    for (int i = 0; i < count_nodes; i++) free(NE[i]);
    free(NE);
}


void graph2matrix(int K, int offset_elements[], int elements[], int count_nodes, int array[]) {
    int left, right, row, col;

    for (int k = 0; k < K; k++) {
        int i = offset_elements[k];
        for (; i < offset_elements[k+1] - 1; ++i) {
            left = elements[i];
            right = elements[i+1];

            row = left;
            col = right;

            array[row*count_nodes + col] = 1;

            row = right;
            col = left;

            array[row*count_nodes + col] = 1;

        }

        left = elements[i];
        right = elements[offset_elements[k]];

        row = left;
        col = right;

        array[row*count_nodes + col] = 1;

        row = right;
        col = left;

        array[row*count_nodes + col] = 1;
    }

    for (int i = 0; i < count_nodes; i++) array[i*count_nodes + i] = 1;
}


void csr2matrix(int N, int IA[], int M, int JA[], int array[]) {
    for (int i=0; i < N-1; i++)
        for (int j = IA[i]; j < IA[i+1]; j++)
            array[i*(N-1) + JA[j]] = 1;
}


void matrix2csr(int K, int array[], int N, int IA[], int M, int JA[]) {
    IA[0] = 0;

    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            if (array[i*K + j] == 1) {
                JA[IA[i+1]] = j;
                ++IA[i+1];
            }
}


void solve(int N, int IA[], double x[], int JA[], double A[], double b[], int *k, double *res, int maxiter, double tol, struct time_by_kernel *time) {
    /*
    N - length of IA, x, b
    n - count of iter
    res - norm L2
    */
    double start_time;

    if (DEBUG) printf("\nDEBUG solve:\n");

    double alpha, betta, rho, last_rho = 1;
    double *r = strict_calloc(N, sizeof(double));
    double *z = strict_calloc(N, sizeof(double));
    double *p = strict_malloc(N * sizeof(double));
    double *q = strict_malloc(N * sizeof(double));
    double *M = strict_malloc(N * sizeof(double));

    for (int i = 0; i < N; i++)
        for (int j = IA[i]; j < IA[i+1]; j++)
            if (i == JA[j])
                M[i] = 1/A[j];

    start_time = omp_get_wtime();
    SpMV(N, IA, JA, A, x, r);
    time->SpMV += omp_get_wtime() - start_time;

    start_time = omp_get_wtime();
    axpby(N, -1, r, 1, b);
    time->axpby += omp_get_wtime() - start_time;

    *k = 1;
    while (1) {
        start_time = omp_get_wtime();
        VVbe(N, M, r, z);
        time->VVbe += omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        rho = dot(N, r, z);
        time->dot += omp_get_wtime() - start_time;

        if (*k == 1) {
            copy(N, z, p);
        } else {
            betta = rho / last_rho;

            start_time = omp_get_wtime();
            axpby(N, betta, p, 1, z);
            time->axpby += omp_get_wtime() - start_time;
        }

        start_time = omp_get_wtime();
        SpMV(N, IA, JA, A, p, q);
        time->SpMV += omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        alpha = rho / dot(N, p, q);
        time->dot += omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        axpby(N, 1, x, alpha, p);
        time->axpby += omp_get_wtime() - start_time;

        start_time = omp_get_wtime();
        axpby(N, 1, r, -alpha, q);
        time->axpby += omp_get_wtime() - start_time;

        printf("\nNumber of iteration: %d\tres = %.2f\n", *k, L2(N, r));

        if ((rho < tol*L2(N, b)) || (*k >= maxiter))
            break;
        else
            ++(*k);

        last_rho = rho;
    }

    *res = L2(N, r);

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
    struct time_by_kernel time;

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

    offset_elements = strict_malloc(sizeof(int) * (K+1));
    exepted_memory += sizeof(int) * (K+1);

    elements = strict_malloc(sizeof(int) * K_full_size);
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

        int *array = strict_calloc(count_nodes * count_nodes, sizeof(int));

        graph2matrix(K, offset_elements, elements, count_nodes, array);

        fprintf(f, "\nGraph2Matrix:\n");
        for (int i = 0; i < count_nodes; i++) {
            for (int j = 0; j < count_nodes; j++)
                fprintf(f, "%d ", array[i*count_nodes + j]);
            fprintf(f, "\n");
        }
    }

    int *IA = strict_malloc((count_nodes + 1) * sizeof(int));
    exepted_memory += (count_nodes + 1) * sizeof(int);
    int *JA = strict_malloc((2 * count_edges + count_nodes) * sizeof(int));
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

        int *array = strict_calloc(count_nodes * count_nodes, sizeof(int));

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
    printf("Time of building: %.2f seconds\n", end_time - start_time);
    printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

    free(elements);
    exepted_memory -= sizeof(int) * K_full_size;
    free(offset_elements);
    exepted_memory -= sizeof(int) * (K+1);

    start_time = omp_get_wtime();

    double *A = strict_malloc((2 * count_edges + count_nodes) * sizeof(double));
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

    double *b = strict_malloc(count_nodes * sizeof(double));
    exepted_memory += count_nodes * sizeof(double);

    for (int i = 0; i < count_nodes; i++) b[i] = sin(i);

    end_time = omp_get_wtime();

    printf("\nCalculated values in matrix\n");
    printf("Time of calculating: %.2f seconds\n", end_time - start_time);
    printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

    start_time = omp_get_wtime();

    double *x = strict_calloc(count_nodes, sizeof(double));
    exepted_memory += count_nodes * sizeof(double);

    omp_set_num_threads(count_threads);
    solve(count_nodes, IA, x, JA, A, b, &count_steps, &res, maxiter, tol, &time);

    if (DEBUG) {
        fprintf(f, "\n\nResults:\n");
        fprintf(f, "  res: %.3f\n", res);
        fprintf(f, "  count steps: %d\n", count_steps);
        fprintf(f, "  X: ");
        for (int i = 0; i < count_nodes; i++) fprintf(f, " %.4f", x[i]);
    }

    end_time = omp_get_wtime();

    printf("\nSolved\n");
    printf("Time of calculating: %.2f seconds\n", end_time - start_time);
    printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

    printf("\nTiming by kernels\n");
    printf("Time of dot: %.3f seconds\n", time.dot);
    printf("Time of axpby: %.3f seconds\n", time.axpby);
    printf("Time of SpMV: %.3f seconds\n", time.SpMV);
    printf("Time of VVbe: %.3f seconds\n\n", time.VVbe);

    return 0;
}
