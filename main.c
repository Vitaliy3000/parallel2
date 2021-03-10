#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE_OF_STRING 40
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))


int DEBUG = 0;

char* MAX_SIZE_OF_STRING_ERROR = "Max size of argumet %s must be less 40 letters\n";
char* DUPLICATE_ARGUMENT_ERROR = "Duplicate argument of %s\n";
char* NOT_FOUND_ARGUMENT_ERROR = "Not found argument of %s\n";
const int LOW = -1;
const int BIG = 100000;
char* LOW_ARGUMENT_ERROR = "Argument %s must be bigger -1\n";
char* BIG_ARGUMENT_ERROR = "Argument %s must be lower 100000\n";
char* INTEGER_ARGUMENT_ERROR = "Argument %s must be ineger\n";

char* HELP_STRING = "\
  usage:\n\
    a.exe [--help]\n\
    a.exe --Nx <Nx> --Ny <Ny> --K1 <K1> --K2 <K2> [--file  <filename>]\n\
  options:\n\
    [--help]     Show this screen.\n\
    --Nx         Count of elements in horizontal\n\
    --Ny         Count of elements in vertical\n\
    --K1         Count of squares\n\
    --K2         Count of cut squares\n\
    [--filename] Debug mode\n\
";

int ERROR_FLAG = 0;


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


void parse_int_arg(int argc, char *argv[], char expected_key[], int ptr[]) {
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

    if (*ptr <= LOW) {
        fprintf(stderr, LOW_ARGUMENT_ERROR, expected_key);
        ERROR_FLAG = 1;
        return;
    }

    if (*ptr >= BIG) {
        fprintf(stderr, BIG_ARGUMENT_ERROR, expected_key);
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


void parse_args(int argc, char *argv[], int *Nx, int *Ny, int *K1, int *K2, char *filename) {
    if ((argc == 1) || (strcmp(argv[1], "--help"))==0) {
        printf(HELP_STRING);
        exit(0);
    }

    parse_int_arg(argc, argv, "--Nx", Nx);
    parse_int_arg(argc, argv, "--Ny", Ny);
    parse_int_arg(argc, argv, "--K1", K1);
    parse_int_arg(argc, argv, "--K2", K2);

    parse_str_arg(argc, argv, "--filename", filename);
}


void calc_count_figures(int Nx, int Ny, int K1, int K2, int *count_triangles, int *count_squares) {
    *count_triangles = 2 * ( K2 * ( Nx * Ny / (K1 + K2) ) + MAX(0, Nx * Ny % (K1 + K2) - K1) );
    *count_squares = K1 * ( Nx * Ny / (K1 + K2) ) + MIN(K1, Nx * Ny % (K1 + K2));
}


void gen_graph(int K, int **elements, int Nx, int Ny, int K1, int K2) {
    int SWITCH_FLAG = 0, i = 0, j = 0;

    for (int k = 0; k < K; k++) {
        if (SWITCH_FLAG >= K1) {
            elements[k][0] = i*(Nx+1) + j;
            elements[k][1] = i*(Nx+1) + j + 1;
            elements[k][2] = (i+1)*(Nx+1) + j;
            elements[k][3] = -1;
            elements[k][4] = -1;

            k++;

            elements[k][0] = i*(Nx+1) + j + 1;
            elements[k][1] = (i+1)*(Nx+1) + j + 1;
            elements[k][2] = (i+1)*(Nx+1) + j;
            elements[k][3] = -1;
            elements[k][4] = -1;
        } else {
            elements[k][0] = i*(Nx+1) + j;
            elements[k][1] = i*(Nx+1) + j + 1;
            elements[k][2] = (i+1)*(Nx+1) + j + 1;
            elements[k][3] = (i+1)*(Nx+1) + j;
            elements[k][4] = -1;
        }

        j++;

        if (j == Nx) {
            i++;
            j = 0;
        }

        SWITCH_FLAG++;

        if (SWITCH_FLAG == (K1 + K2)) SWITCH_FLAG = 0;
    }
}


int cmpfunc(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}


void graph2csr(int K, int *elements[], int count_edges, int count_nodes, int *IA,  int *JA) {
    int next, i, j, k, left, right, current_filling, offset;

    int *counter_IA_d = strict_calloc(count_nodes, sizeof(int));

    for (k = 0; k < K; k++) {
        for (i = 0; elements[k][i+1] != -1; i++) {
            left = elements[k][i];
            right = elements[k][i+1];

            ++counter_IA_d[left];
            ++counter_IA_d[right];
        }

        left = elements[k][i];
        right = elements[k][0];
    
        ++counter_IA_d[left];
        ++counter_IA_d[right];
    }

    int **NE_d = strict_malloc(count_nodes * sizeof(int*));
    for (i = 0; i < count_nodes; i++) NE_d[i] = strict_malloc(counter_IA_d[i] * sizeof(int));
    for (i = 0; i < count_nodes; i++) counter_IA_d[i] = 0;

    for (k = 0; k < K; k++) {
        for (i = 0; elements[k][i+1] != -1; i++) {
            left = elements[k][i];
            right = elements[k][i+1];

            NE_d[left][counter_IA_d[left]] = right;
            NE_d[right][counter_IA_d[right]] = left;

            ++counter_IA_d[left];
            ++counter_IA_d[right];
        }

        left = elements[k][i];
        right = elements[k][0];

        NE_d[left][counter_IA_d[left]] = right;
        NE_d[right][counter_IA_d[right]] = left;

        ++counter_IA_d[left];
        ++counter_IA_d[right];
    }


    for (i = 0; i < count_nodes; i++) qsort(NE_d[i], counter_IA_d[i], sizeof(int), cmpfunc);

    int *counter_IA = strict_calloc(count_nodes, sizeof(int));

    for (i = 0; i < count_nodes; i++) {
        for (j = 0; j < counter_IA_d[i]-1; j++) {
            if (NE_d[i][j] != NE_d[i][j+1]) {
                ++counter_IA[i];
            }
        }

        ++counter_IA[i];
    }

    int **NE = strict_malloc(count_nodes * sizeof(int*));
    for (i = 0; i < count_nodes; i++) NE[i] = strict_malloc(counter_IA[i] * sizeof(int));
    for (i = 0; i < count_nodes; i++) counter_IA[i] = 0;

    for (i = 0; i < count_nodes; i++) {
        for (j = 0; j < counter_IA_d[i]-1; j++) {
            if (NE_d[i][j] != NE_d[i][j+1]) {
                NE[i][counter_IA[i]] = NE_d[i][j];
                ++counter_IA[i];
            }
        }

        NE[i][counter_IA[i]] = NE_d[i][j];
        ++counter_IA[i];
    }

    IA[0] = 0;
    for (i = 1; i < count_nodes+1; i++) IA[i] = IA[i-1] + counter_IA[i-1];
    for (i = 0; i < count_nodes; i++) counter_IA[i] = 0;

    k = 0;
    for (i = 0; i < count_nodes; i++) {
        for (j = 0; j < IA[i+1] - IA[i]; j++) {
            JA[k] = NE[i][j];
            ++k;
        }
    }
}


void graph2matrix(int K, int *elements[], int count_nodes, int *array) {
    int i, k, left, right, row, col;

    for (k = 0; k < K; k++) {
        for (i = 0; elements[k][i+1] != -1; i++) {
            left = elements[k][i];
            right = elements[k][i+1];

            row = left;
            col = right;

            array[row*count_nodes + col] = 1;

            row = right;
            col = left;

            array[row*count_nodes + col] = 1;
        }

        left = elements[k][i];
        right = elements[k][0];

        row = left;
        col = right;

        array[row*count_nodes + col] = 1;

        row = right;
        col = left;

        array[row*count_nodes + col] = 1;
    }
}


void csr2matrix(int N, int IA[], int M, int JA[], int array[]) {
    int i, j;

    for (i=0; i < N-1; i++)
        for (j = IA[i]; j < IA[i+1]; j++)
            array[i*(N-1) + JA[j]] = 1;
}


void matrix2csr(int K, int array[], int N, int IA[], int M, int JA[]) {
    int i, j;

    IA[0] = 0;

    for (i = 0; i < K; i++)
        for (j = 0; j < K; j++)
            if (array[i*K + j] == 1) {
                JA[IA[i+1]] = j;
                ++IA[i+1];
            }
}


int main(int argc, char *argv[]) {
    int i, j, k, Nx, Ny, K1, K2, K, count_triangles, count_squares;
    int **elements;
    FILE *f;
    char filename[MAX_SIZE_OF_STRING] = "";

    parse_args(argc, argv, &Nx, &Ny, &K1, &K2, filename);

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
        fprintf(f, "file = %s\n", filename);
    }

    calc_count_figures(Nx, Ny, K1, K2, &count_triangles, &count_squares);
    K = count_squares + count_triangles;

    elements = strict_malloc(sizeof(int*) * K);

    for (i = 0; i < K; i++)
        elements[i] = strict_malloc(sizeof(int) * 5);

    gen_graph(K, elements, Nx, Ny, K1, K2);

    int count_edges = Ny*(Nx+1) + Nx*(Ny+1) + count_triangles / 2;
    int count_nodes = (Nx+1) * (Ny+1);

    if (DEBUG) {
        fprintf(f, "\nGraph:\n");
        for (i = 0; i < K; i++) {
            fprintf(f, "    Element:");
            for (j = 0; j < 4; j++)
                fprintf(f, " %d", elements[i][j]);
            fprintf(f, "\n");
        }

        fprintf(f, "\n");

        fprintf(f, "    Count nodes: %d\n", count_nodes);
        fprintf(f, "    Count edges: %d\n", count_edges);

        int *array = strict_calloc(count_nodes * count_nodes, sizeof(int));

        graph2matrix(K, elements, count_nodes, array);

        fprintf(f, "\nGraph2Matrix:\n");
        for (i = 0; i < count_nodes; i++) {
            for (j = 0; j < count_nodes; j++)
                fprintf(f, "%d ", array[i*count_nodes + j]);
            fprintf(f, "\n");
        }
    }

    int *IA = strict_malloc((count_nodes + 1) * sizeof(int));
    int *JA = strict_malloc((2 * count_edges) * sizeof(int));

    graph2csr(K, elements, count_edges, count_nodes, IA, JA);

    if (DEBUG) {
        fprintf(f, "\nCSR:");
        fprintf(f, "\n    IA:");
        for (i = 0; i < count_nodes + 1; i++)
            fprintf(f, " %d", IA[i]);
        fprintf(f, "\n    JA:");
        for (i = 0; i < 2 * count_edges; i++)
            fprintf(f, " %d", JA[i]);

        int *array = strict_calloc(count_nodes * count_nodes, sizeof(int));

        csr2matrix(count_nodes + 1, IA, 2 * count_edges, JA, array);

        fprintf(f, "\n\nCSR2Matrix:\n");
        for (i = 0; i < count_nodes; i++) {
            for (j = 0; j < count_nodes; j++)
                fprintf(f, "%d ", array[i*count_nodes + j]);
            fprintf(f, "\n");
        }
    }

    return 0;
}
