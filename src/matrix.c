#include "matrix.h"


int cmpfunc(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}


void graph2csr(int K, int offset_elements[], int elements[], int count_edges, int count_nodes, int IA[],  int JA[]) {
    int left, right, flag;

    int *counter_IA_d = (int*) strict_calloc(count_nodes, sizeof(int));

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

    int **NE_d = (int**) strict_malloc(count_nodes * sizeof(int*));
    for (int i = 0; i < count_nodes; i++) NE_d[i] = (int*) strict_malloc(counter_IA_d[i] * sizeof(int));
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

    int *counter_IA = (int*) strict_calloc(count_nodes, sizeof(int));

    for (int i = 0; i < count_nodes; i++) {
        for (int j = 0; j < counter_IA_d[i]-1; j++) {
            if (NE_d[i][j] != NE_d[i][j+1]) {
                ++counter_IA[i];
            }
        }

        ++counter_IA[i];
    }

    int **NE = (int**) strict_malloc(count_nodes * sizeof(int*));
    for (int i = 0; i < count_nodes; i++) NE[i] = (int*) strict_malloc(counter_IA[i] * sizeof(int));
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
