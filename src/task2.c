#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
// #include <omp.h>
#include "parser.h"
#include "lib.h"
#include "matrix.h"
#include "kernals.h"


#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))


// extern int ERROR_FLAG;

int DEBUG = 0;
int ROOT = 0;

int size, rank;

int *L2G;

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
const int LOW_Px = 1;
const int HIGH_Px = 500;
const int LOW_Py = 1;
const int HIGH_Py = 500;

char* HELP_STRING = "\
  usage:\n\
    mpirun -np <np> task2 [--help]\n\
    mpirun -np <np> task2 --Nx <Nx> --Ny <Ny> --K1 <K1> --K2 <K2> --Px <Px> --Py <Py> --maxiter <N> --tol <tol> [--file  <filename>]\n\
  options:\n\
    [--help]     Show this screen.\n\
    --Nx         Count of elements in horizontal\n\
    --Ny         Count of elements in vertical\n\
    --K1         Count of squares\n\
    --K2         Count of cut squares\n\
    --Px         Count of parts by x\n\
    --Py         Count of parts by y\n\
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


void mpi_exit(int arg) {
    MPI_Finalize();
    exit(arg);
}


void parse_args(int argc, char *argv[], int *Nx, int *Ny, int *K1, int *K2, int *Px, int *Py, int *maxiter, double *tol, char *filename) {
    parse_int_arg(argc, argv, "--Nx", Nx, LOW_Nx, HIGH_Nx);
    parse_int_arg(argc, argv, "--Ny", Ny, LOW_Ny, HIGH_Ny);
    parse_int_arg(argc, argv, "--K1", K1, LOW_K1, HIGH_K1);
    parse_int_arg(argc, argv, "--K2", K2, LOW_K2, HIGH_K2);
    parse_int_arg(argc, argv, "--Px", Px, LOW_Px, HIGH_Px);
    parse_int_arg(argc, argv, "--Py", Py, LOW_Py, HIGH_Py);
    parse_int_arg(argc, argv, "--maxiter", maxiter, LOW_maxiter, HIGH_maxiter);
    parse_double_arg(argc, argv, "--tol", tol, LOW_tol, HIGH_tol);

    parse_str_arg(argc, argv, "--filename", filename);
}


void calc_count_figures(int Nx, int Ny, int K1, int K2, int ib, int ie, int jb, int je, int *count_triangles, int *count_squares) {
    *count_squares = 0;
    *count_triangles = 0;

    for (int i = ib; i < ie-1; i++) {
        for (int j = jb; j < je-1; j++) {
            if ((i*Nx + j) % (K1+K2) < K1)
                (*count_squares)++;
            else 
                (*count_triangles) += 2;
        }
    }
}


int get_local_number(int i, int j, int ib, int ie, int jb, int je) {
    return (i - ib) * (je - jb) + (j - jb);
}

void gen_graph(int offset_elements[], int elements[], int Nx, int Ny, int K1, int K2, int ib, int ie, int jb, int je) {
    int k=0, p=0;

    offset_elements[0] = 0;

    for (int i = ib; i < ie-1; i++) {
        for (int j = jb; j < je-1; j++) {
            if ((i*Nx + j) % (K1+K2) < K1) {
                elements[p++] = get_local_number(i, j, ib, ie, jb, je);
                elements[p++] = get_local_number(i, j+1, ib, ie, jb, je);
                elements[p++] = get_local_number(i+1, j+1, ib, ie, jb, je);
                elements[p++] = get_local_number(i+1, j, ib, ie, jb, je);
                offset_elements[k+1] = offset_elements[k] + 4;
                ++k;
            } else {
                elements[p++] = get_local_number(i, j, ib, ie, jb, je);
                elements[p++] = get_local_number(i, j+1, ib, ie, jb, je);
                elements[p++] = get_local_number(i+1, j, ib, ie, jb, je);

                offset_elements[k+1] = offset_elements[k] + 3;
                ++k;

                elements[p++] = get_local_number(i, j+1, ib, ie, jb, je);
                elements[p++] = get_local_number(i+1, j+1, ib, ie, jb, je);
                elements[p++] = get_local_number(i+1, j, ib, ie, jb, je);

                offset_elements[k+1] = offset_elements[k] + 3;
                ++k;
            }
        }
    }
}

int get_global_number(int i, int j, int Nx) {
    return i * (Nx+1) + j;
}

int get_number_proc(int i, int j, int Nx, int Ny, int Px, int Py) {
    int on_proc_x = (Nx+1) / Px;
    int over_proc_x = (Nx+1) % Px;
    int on_proc_y = (Ny+1) / Py;
    int over_proc_y = (Ny+1) % Py;

    if (j > (on_proc_x + 1) * over_proc_x )
        j = ( j - (on_proc_x + 1) * over_proc_x ) / on_proc_x + over_proc_x;
    else
        j = j / (on_proc_x + 1);

    if (i > (on_proc_y + 1) * over_proc_y )
        i = ( i - (on_proc_y + 1) * over_proc_y ) / on_proc_y + over_proc_y;
    else
        i = i / (on_proc_y + 1);

    return i * Px + j;
}


void update_galo(int N, double x[], int offset_on_proc[], int indexes_on_send[], int indexes_on_recv[]) {
    int count_neighbours = 0;
    for (int i = 0; i < size; i++) count_neighbours += offset_on_proc[i+1] != offset_on_proc[i];

    MPI_Status *statuses = (MPI_Status*) strict_malloc(sizeof(MPI_Status) * 2 * count_neighbours);
    MPI_Request *requests = (MPI_Request*) strict_malloc(sizeof(MPI_Request) * 2 * count_neighbours);

    double *recv = (double*) strict_malloc(sizeof(double) * offset_on_proc[size]);
    double *send = (double*) strict_malloc(sizeof(double) * offset_on_proc[size]);

    int ncom = 0;

    for (int i = 0; i < size; i++) {
        int count_on_recv = offset_on_proc[i+1] - offset_on_proc[i];
        if (count_on_recv > 0) {
            MPI_Irecv(&recv[offset_on_proc[i]], count_on_recv, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &requests[ncom]);
            ++ncom;
        }
    }

    for (int i = 0; i < size; i++)
        for (int j = offset_on_proc[i]; j < offset_on_proc[i+1]; j++) {
            send[j] = x[indexes_on_send[j]];
        }


    for (int i = 0; i < size; i++) {
        int count_on_send = offset_on_proc[i+1] - offset_on_proc[i];
        if (count_on_send > 0) {
            MPI_Isend(&send[offset_on_proc[i]], count_on_send, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &requests[ncom]);
            ++ncom;
        }
    }

    MPI_Waitall(count_neighbours, requests, statuses);

    for (int i = 0; i < size; i++)
        for (int j = offset_on_proc[i]; j < offset_on_proc[i+1]; j++) {
            x[indexes_on_recv[j]] = recv[j];
        }
}

double dot_MPI(int N, double x[], double y[]) {
    double res = dot(N, x, y);
    MPI_Allreduce(&res, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return res;
}

double L2_MPI(int N, double x[]) {
    return sqrtf(dot_MPI(N, x, x));
}

void solve(int N, int No, int IA[], double x[], int JA[], double A[], double b[], int *k, double *res, int maxiter, double tol, int offset_on_proc[], int indexes_on_send[], int indexes_on_recv[], struct res_by_kernel *gflops) {
    /*
    N - length of IA - 1, x, b
    No - length own nodes
    k - count of iter
    res - norm L2
    */
    // double start_time;

    gflops->gflops_dot = 0;
    gflops->gflops_axpby = 0;
    gflops->gflops_SpMV = 0;
    gflops->gflops_VVbe = 0;

    gflops->time_dot = 0;
    gflops->time_axpby = 0;
    gflops->time_SpMV = 0;
    gflops->time_VVbe = 0;

    if ((DEBUG) && (rank == ROOT)) printf("\nDEBUG solve:\n");

    double alpha, betta, rho, last_rho = 1, start_time;
    double *r = (double*) strict_calloc(N, sizeof(double));
    double *z = (double*) strict_calloc(No, sizeof(double));
    double *p = (double*) strict_malloc(N * sizeof(double));
    double *q = (double*) strict_malloc(N * sizeof(double));
    double *M = (double*) strict_malloc(No * sizeof(double));

    for (int i = 0; i < No; i++)
        for (int j = IA[i]; j < IA[i+1]; j++)
            if (i == JA[j])
                M[i] = 1/A[j];

    start_time = MPI_Wtime();
    update_galo(N, x, offset_on_proc, indexes_on_send, indexes_on_recv);
    SpMV(N, IA, JA, A, x, r);
    gflops->time_SpMV += MPI_Wtime() - start_time;
    gflops->gflops_SpMV += 2 * 1E-9 * IA[N];

    start_time = MPI_Wtime();
    axpby(No, -1, r, 1, b);
    gflops->time_axpby += MPI_Wtime() - start_time;
    gflops->gflops_axpby += 3 * 1E-9 * N;

    *k = 1;
    while (1) {
        start_time = MPI_Wtime();
        VVbe(No, M, r, z);
        gflops->time_VVbe += MPI_Wtime() - start_time;
        gflops->gflops_VVbe += 1E-9 * N;

        start_time = MPI_Wtime();
        rho = dot_MPI(No, r, z);
        gflops->time_dot += MPI_Wtime() - start_time;
        gflops->gflops_dot += 2 * 1E-9 * N;

        if (*k == 1) {
            copy(No, z, p);
        } else {
            betta = rho / last_rho;

            start_time = MPI_Wtime();
            axpby(No, betta, p, 1, z);
            gflops->time_axpby += MPI_Wtime() - start_time;
            gflops->gflops_axpby += 3 * 1E-9 * N;
        }

        start_time = MPI_Wtime();
        update_galo(N, p, offset_on_proc, indexes_on_send, indexes_on_recv);
        SpMV(N, IA, JA, A, p, q);
        gflops->time_SpMV += MPI_Wtime() - start_time;
        gflops->gflops_SpMV += 2 * 1E-9 * IA[N];

        start_time = MPI_Wtime();
        alpha = rho / dot_MPI(No, p, q);
        *res = L2_MPI(No, q);

        gflops->time_dot += MPI_Wtime() - start_time;
        gflops->gflops_dot += 2 * 1E-9 * N;

        start_time = MPI_Wtime();
        axpby(No, 1, x, alpha, p);
        gflops->time_axpby += MPI_Wtime() - start_time;
        gflops->gflops_axpby += 3 * 1E-9 * N;

        start_time = MPI_Wtime();
        axpby(No, 1, r, -alpha, q);
        gflops->time_axpby += MPI_Wtime() - start_time;
        gflops->gflops_axpby += 3 * 1E-9 * N;

        *res = L2_MPI(No, r);

        if (rank == ROOT)
            printf("\nNumber of iteration: %d\tres = %.5f\n", *k, *res);

        if ((*res / L2_MPI(No, b) < tol) || (*k >= maxiter))
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


void iswap(int *a, int *b) {
    (*a) = (*a) + (*b);
    (*b) = (*a) - (*b);
    (*a) = (*a) - (*b);
}


int cmpfunc_l2g(const void * a, const void * b) {
    return ( L2G[*(int*)a] - L2G[*(int*)b] );
}


int main(int argc, char *argv[]) {
    int exepted_memory = 0;
    double start_time, end_time;

    int errCode;

    if ((errCode = MPI_Init(&argc, &argv)) != 0) {
        perror("MPI_Init fail");
        return errCode;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if ((argc == 1) || (strcmp(argv[1], "--help"))==0) {
        if (rank == ROOT) printf(HELP_STRING);
        mpi_exit(0);
    }

    int Nx, Ny, Px, Py, K1, K2, maxiter;
    double tol;
    FILE *f;
    char filename[MAX_SIZE_OF_STRING] = "";

    parse_args(argc, argv, &Nx, &Ny, &K1, &K2, &Px, &Py, &maxiter, &tol, filename);

    if (ERROR_FLAG) mpi_exit(EXIT_FAILURE);
    if (size != (Px*Py)) {
        perror("count procs must be equal Px * Py");
        mpi_exit(EXIT_FAILURE);
    }

    if (filename[0] != '\0') {
        f = fopen(filename, "w");
        DEBUG = 1;
    }

    start_time = MPI_Wtime();

    int cor_x = rank % Px;
    int cor_y = rank / Px;

    int count_galo_nodes_x = 0;
    int count_galo_nodes_y = 0;

    int jb = ((Nx+1) / Px) * cor_x + MIN((Nx+1) % Px, cor_x);
    int je = jb + ((Nx+1) / Px) + (cor_x < (Nx+1) % Px);
    count_galo_nodes_x += (jb != 0);
    jb -= (jb != 0);
    count_galo_nodes_x += (je != (Nx+1));
    je += (je != (Nx+1));

    int ib = ((Ny+1) / Py) * cor_y + MIN((Ny+1) % Py, cor_y);
    int ie = ib + ((Ny+1) / Py) + (cor_y < (Ny+1) % Py);
    count_galo_nodes_y += (ib!=0);
    ib -= (ib!=0);
    count_galo_nodes_y += (ie != (Ny+1));
    ie += (ie != (Ny+1));

    int count_galo_nodes = count_galo_nodes_x * (ie - ib) + count_galo_nodes_y * (je - jb) - count_galo_nodes_x*count_galo_nodes_y;
    int N = (ie-ib) * (je-jb);
    int No = N - count_galo_nodes;

    if (No == 0) {
        N = 0;
        ib = 0;
        ie = 0;
        jb = 0;
        je = 0;
    }

    L2G = (int*) strict_realloc(L2G, sizeof(int) * N);
    int *Part = (int*) strict_malloc(sizeof(int) * N);
    int *SORTED_L2G = (int*) strict_malloc(sizeof(int) * N);

    for (int i = ib; i < ie; i++)
        for (int j = jb; j < je; j++) {
            L2G[get_local_number(i, j, ib, ie, jb, je)] = get_global_number(i, j, Nx);
            Part[get_local_number(i, j, ib, ie, jb, je)] = get_number_proc(i, j, Nx, Ny, Px, Py);
            SORTED_L2G[get_local_number(i, j, ib, ie, jb, je)] = get_local_number(i, j, ib, ie, jb, je);
        }

    int rb = No;
    for (int lb = 0; lb < No; lb++) {
        if (Part[lb] != rank) {
            while (Part[rb] != rank) ++rb;
            SORTED_L2G[lb] = rb;
            SORTED_L2G[rb] = lb;
            ++rb;
        }
    }

    if (DEBUG) {
        if (rank == ROOT) {
            fprintf(f, "Grid:\n");
            for (int i = 0; i < (Ny+1); i++) {
                fprintf(f, "  ");
                for (int j = 0; j < (Nx+1); j++)
                    fprintf(f, "%d ", i*(Nx+1) + j);
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
    }

    if (DEBUG) {
        int nlx = (je-jb);
        int nly = (ie-ib);
        if (rank == ROOT) {
            int *ngx = (int*) strict_malloc(sizeof(int) * size);
            int *ngy = (int*) strict_malloc(sizeof(int) * size);
            MPI_Gather(&nlx, 1, MPI_INT, ngx, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&nly, 1, MPI_INT, ngy, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            int *recvcounts = strict_malloc(sizeof(int) * size);
            for (int i = 0; i < size; i++) recvcounts[i] = ngx[i] * ngy[i];
            int *displs = (int*) strict_malloc(sizeof(int) * size);
            displs[0] = 0;
            for (int i = 1; i < size; i++) displs[i] = displs[i-1] + recvcounts[i-1];
            int *L2G_ALL = (int*) strict_malloc(sizeof(int) * (displs[size-1] + recvcounts[size-1]));
            int *Part_ALL = (int*) strict_malloc(sizeof(int) * (displs[size-1] + recvcounts[size-1]));
            int *SORTED_L2G_ALL = (int*) strict_malloc(sizeof(int) * (displs[size-1] + recvcounts[size-1]));
            MPI_Gatherv(L2G, N, MPI_INT, L2G_ALL, recvcounts, displs, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(Part, N, MPI_INT, Part_ALL, recvcounts, displs, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(SORTED_L2G, N, MPI_INT, SORTED_L2G_ALL, recvcounts, displs, MPI_INT, ROOT, MPI_COMM_WORLD);

            fprintf(f, "Local & Global & Part grids:\n");
            for (int i = 0; i < size; i++) {
                fprintf(f, "rank = %d: (%d, %d): \n", i, ngx[i], ngy[i]);
                for (int k = 0; k < ngy[i]; k++) {
                    for (int j = 0; j < ngx[i]; j++) fprintf(f, "%d ", k*ngx[i] + j);
                    fprintf(f, "\t");
                    for (int j = 0; j < ngx[i]; j++) fprintf(f, "%d ", L2G_ALL[displs[i] + k*ngx[i] + j]);
                    fprintf(f, "\t");
                    for (int j = 0; j < ngx[i]; j++) fprintf(f, "%d ", Part_ALL[displs[i] + k*ngx[i] + j]);
                    fprintf(f, "\n");
                }
                fprintf(f, "\n");
            }

            fprintf(f, "SORTED_L2G_ALL:\n");
            for (int i = 0; i < size; i++) {
                fprintf(f, "rank = %d: (%d, %d): \n", i, ngx[i], ngy[i]);
                for (int k = 0; k < ngy[i]; k++) {
                    for (int j = 0; j < ngx[i]; j++) fprintf(f, "%d->%d ", k*ngx[i] + j, SORTED_L2G_ALL[displs[i] + k*ngx[i] + j]);
                    fprintf(f, "\n");
                }
                fprintf(f, "\n");
            }


            fprintf(f, "\n");
        } else {
            MPI_Gather(&nlx, 1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&nly, 1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(L2G, N, MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(Part, N, MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(SORTED_L2G, N, MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
        }
    }

    int count_triangles, count_squares;
    calc_count_figures(Nx, Ny, K1, K2, ib, ie, jb, je, &count_triangles, &count_squares);

    int *offset_elements = (int*) strict_malloc(sizeof(int) * (count_triangles + count_squares + 1));
    exepted_memory += sizeof(int) * (count_triangles + count_squares + 1);

    int *elements = (int*) strict_malloc(sizeof(int) * (3*count_triangles + 4*count_squares));
    exepted_memory += sizeof(int) * (3*count_triangles + 4*count_squares);

    gen_graph(offset_elements, elements, Nx, Ny, K1, K2, ib, ie, jb, je);

    if (DEBUG) {
        int N_offset = count_triangles + count_squares + 1;
        int N_elements = 3*count_triangles + 4*count_squares;

        if (rank == ROOT) {
            int *sizes_of_offsets = (int*) strict_malloc(sizeof(int) * size);
            int *sizes_of_elements = (int*) strict_malloc(sizeof(int) * size);

            MPI_Gather(&N_offset, 1, MPI_INT, sizes_of_offsets, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&N_elements, 1, MPI_INT, sizes_of_elements, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

            int N_offset_all = 0;
            for (int i = 0; i < size; i++) N_offset_all += sizes_of_offsets[i];

            int N_elements_all = 0;
            for (int i = 0; i < size; i++) N_elements_all += sizes_of_elements[i];

            int *all_offsets = (int*) strict_malloc(sizeof(int*) * N_offset_all);
            int *all_elements = (int*) strict_malloc(sizeof(int*) * N_elements_all);

            int *recvcounts_offsets = strict_malloc(sizeof(int) * size);
            for (int i = 0; i < size; i++) recvcounts_offsets[i] = sizes_of_offsets[i];
            int *displs_offsets = strict_malloc(sizeof(int) * size);
            displs_offsets[0] = 0;
            for (int i = 1; i < size; i++) displs_offsets[i] = displs_offsets[i-1] + recvcounts_offsets[i-1];

            int *recvcounts_elements = strict_malloc(sizeof(int) * size);
            for (int i = 0; i < size; i++) recvcounts_elements[i] = sizes_of_elements[i];
            int *displs_elements = strict_malloc(sizeof(int) * size);
            displs_elements[0] = 0;
            for (int i = 1; i < size; i++) displs_elements[i] = displs_elements[i-1] + recvcounts_elements[i-1];

            MPI_Gatherv(offset_elements, N_offset, MPI_INT, all_offsets, recvcounts_offsets, displs_offsets, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(elements, N_elements, MPI_INT, all_elements, recvcounts_elements, displs_elements, MPI_INT, ROOT, MPI_COMM_WORLD);

            for (int r = 0; r < size; r++) {
                fprintf(f, "\nLocal graph (rank=%d):\n", r);
                for (int j = 0; j < sizes_of_offsets[r]-1; ++j){
                    fprintf(f, "    Element(%d): ", j);
                    for (int i = all_offsets[displs_offsets[r]+j]; i < all_offsets[displs_offsets[r]+j+1]; ++i)
                        fprintf(f, " %d", all_elements[displs_elements[r]+i]);

                    fprintf(f, "\n");
                }
            }

            fprintf(f, "\n");
        } else {
            MPI_Gather(&N_offset, 1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&N_elements, 1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);

            MPI_Gatherv(offset_elements, N_offset, MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(elements, N_elements, MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
        }
    }

    for (int i = 0; i < No; i++)
        if (SORTED_L2G[i] != i) {
            iswap(&L2G[i], &L2G[SORTED_L2G[i]]);
            iswap(&Part[i], &Part[SORTED_L2G[i]]);
        }

    for (int i = 0; i < (3*count_triangles + 4*count_squares); i++) elements[i] = SORTED_L2G[elements[i]];

    int count_edges = (je-jb-1)*(ie-ib) + (je-jb)*(ie-ib-1) + count_triangles / 2;
    int count_nodes = N;

    int *IA = (int*) strict_malloc((N + 1) * sizeof(int));
    exepted_memory += (N + 1) * sizeof(int);
    int *JA = (int*) strict_malloc((2 * count_edges + count_nodes) * sizeof(int));
    exepted_memory += (2 * count_edges + count_nodes) * sizeof(int);

    graph2csr(count_triangles + count_squares, offset_elements, elements, count_edges, count_nodes, IA, JA);

    end_time = MPI_Wtime();

    free(elements);
    exepted_memory -= sizeof(int) * (3*count_triangles + 4*count_squares);
    free(offset_elements);
    exepted_memory -= sizeof(int) * (count_triangles + count_squares+1);

    if (rank == ROOT) {
        printf("\nGraph in CSR format was build\n");
        printf("Time of building: %.5f seconds\n", end_time - start_time);
        // printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);
    }

    start_time = MPI_Wtime();

    double *A = (double*) strict_malloc(IA[N] * sizeof(double));
    exepted_memory += IA[N] * sizeof(double);

    int k = 0;
    for (int i = 0; i < N; i++) {
        double sum = 0;
        int flag;
        for (int j = IA[i]; j < IA[i+1]; j++) {
            if (i != JA[j]) {
                A[k] = cos(L2G[i]*L2G[JA[j]] + L2G[i] + L2G[JA[j]]);
                // printf("%d, %f %d\n",k, A[k], L2G[i]*L2G[JA[j]] + L2G[i] + L2G[JA[j]]);
                sum += fabs(A[k]);
            } else {
                flag = k;
            }
            ++k;
        }
        A[flag] = 1.234*sum;
    }

    double *b = (double*) strict_malloc(No * sizeof(double));
    exepted_memory += N * sizeof(double);

    for (int i = 0; i < No; i++) b[i] = sin(L2G[i]);

    double *x = (double*) strict_calloc(N, sizeof(double));
    exepted_memory += N * sizeof(double);

    end_time = MPI_Wtime();

    if (rank == ROOT) {
        printf("\nCalculated values in matrix\n");
        printf("Time of calculating: %.5f seconds\n", end_time - start_time);
        // printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);
    }

    start_time = MPI_Wtime();

    int *count_on_process = (int*) strict_calloc(size, sizeof(int));

    for (int i = 0; i < N; ++i)
        for (int j = IA[i]; j < IA[i+1]; j++)
            if ((Part[i] == rank) && (Part[JA[j]] != rank))
                ++count_on_process[Part[JA[j]]];

    int s = 0;
    for (int i = 0; i < size; i++) s += count_on_process[i];
    int *indexes_on_send = (int*) strict_malloc(s * sizeof(int));
    int *indexes_on_recv = (int*) strict_malloc(s * sizeof(int));
    int *offset_on_proc = (int*) strict_malloc((size+1) * sizeof(int));

    offset_on_proc[0] = 0;

    for (int i = 1; i < (size+1); i++) offset_on_proc[i] = offset_on_proc[i-1] + count_on_process[i-1];

    for (int i = 0; i < size; i++) count_on_process[i] = 0;

    for (int i = 0; i < N; ++i)
        for (int j = IA[i]; j < IA[i+1]; j++)
            if ((Part[i] == rank) && (Part[JA[j]] != rank)) {
                int n_proc = Part[JA[j]];
                indexes_on_send[offset_on_proc[n_proc] + count_on_process[n_proc]] = i;
                indexes_on_recv[offset_on_proc[n_proc] + count_on_process[n_proc]] = JA[j];
                ++count_on_process[n_proc];
            }

    free(count_on_process);


    for (int i = 0; i < size; i++) {
        int count = offset_on_proc[i+1]-offset_on_proc[i];
        if (count > 0) {
            qsort(&indexes_on_send[offset_on_proc[i]], count, sizeof(int), cmpfunc_l2g);
            qsort(&indexes_on_recv[offset_on_proc[i]], count, sizeof(int), cmpfunc_l2g);
        }
    }

    if (rank == ROOT) {
        printf("\nBuilt scheme of communication\n");
        printf("Time of building: %.5f seconds\n", end_time - start_time);
        // printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);
    }

    if (DEBUG) {
        if (rank == ROOT) {
            int *N_ALL = (int*) strict_malloc(sizeof(int) * size);
            int *No_ALL = (int*) strict_malloc(sizeof(int) * size);

            MPI_Gather(&N, 1, MPI_INT, N_ALL, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&No, 1, MPI_INT, No_ALL, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

            int *recvcounts = strict_malloc(sizeof(int) * size);
            for (int i = 0; i < size; i++) recvcounts[i] = N_ALL[i];
            int *displs = (int*) strict_malloc(sizeof(int) * size);
            displs[0] = 0;
            for (int i = 1; i < size; i++) displs[i] = displs[i-1] + recvcounts[i-1];
            int *L2G_ALL = (int*) strict_malloc(sizeof(int) * (displs[size-1] + recvcounts[size-1]));
            MPI_Gatherv(L2G, N, MPI_INT, L2G_ALL, recvcounts, displs, MPI_INT, ROOT, MPI_COMM_WORLD);

            for (int r = 0; r < size; r++) {
                fprintf(f, "\nL2G(rank=%d):", r);
                for (int j = 0; j < No_ALL[r]; ++j)
                    fprintf(f, " %d", L2G_ALL[displs[r]+j]);

                fprintf(f, " | ");

                for (int j = No_ALL[r]; j < N_ALL[r]; ++j)
                    fprintf(f, " %d", L2G_ALL[displs[r]+j]);

                fprintf(f, "\n");
            }

            int *offset_on_proc_all = strict_malloc(sizeof(int) * size * (size+1));

            MPI_Gather(offset_on_proc, size+1, MPI_INT, offset_on_proc_all, size+1, MPI_INT, ROOT, MPI_COMM_WORLD);

            int *recvcounts_indexes = strict_malloc(sizeof(int) * size);
            for (int r = 0; r < size; r++) recvcounts_indexes[r] = offset_on_proc_all[(size+1)*r+size];
            int *displs_indexes = (int*) strict_malloc(sizeof(int) * size);
            displs_indexes[0] = 0;
            for (int r = 1; r < size; r++) displs_indexes[r] = displs_indexes[r-1] + recvcounts_indexes[r-1];

            int *indexes_on_send_all = (int*) strict_malloc(sizeof(int) * (displs_indexes[size-1] + recvcounts_indexes[size-1]));
            int *indexes_on_recv_all = (int*) strict_malloc(sizeof(int) * (displs_indexes[size-1] + recvcounts_indexes[size-1]));
            MPI_Gatherv(indexes_on_send, offset_on_proc[size], MPI_INT, indexes_on_send_all, recvcounts_indexes, displs_indexes, MPI_INT, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(indexes_on_recv, offset_on_proc[size], MPI_INT, indexes_on_recv_all, recvcounts_indexes, displs_indexes, MPI_INT, ROOT, MPI_COMM_WORLD);

            for (int r = 0; r < size; r++) {
                fprintf(f, "\nS&R(rank=%d):\n", r);

                for (int i = 0; i < size; i++) {
                    if (offset_on_proc_all[(size+1)*r + i] != offset_on_proc_all[(size+1)*r + i + 1]) {
                        fprintf(f, " send nums: ");
                        for (int j = offset_on_proc_all[(size+1)*r + i]; j < offset_on_proc_all[(size+1)*r+i+1]; j++)
                            fprintf(f, "%d ", indexes_on_send_all[displs_indexes[r] + j]);
                        fprintf(f, " on rank=%d\n", i);
                    }
                }

                for (int i = 0; i < size; i++) {
                    if (offset_on_proc_all[(size+1)*r + i] != offset_on_proc_all[(size+1)*r + i + 1]) {
                        fprintf(f, " recv nums: ");
                        for (int j = offset_on_proc_all[(size+1)*r + i]; j < offset_on_proc_all[(size+1)*r+i+1]; j++)
                            fprintf(f, "%d ", indexes_on_recv_all[displs_indexes[r] + j]);
                        fprintf(f, " on rank=%d\n", i);
                    }
                }

                fprintf(f, "\n\n");
            }
        } else {
            MPI_Gather(&N, 1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gather(&No, 1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(L2G, N, MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gather(offset_on_proc, size+1, MPI_INT, NULL, 0, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(indexes_on_send, offset_on_proc[size], MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
            MPI_Gatherv(indexes_on_recv, offset_on_proc[size], MPI_INT, NULL, NULL, NULL, NULL, ROOT, MPI_COMM_WORLD);
        }
    }

    int count_steps;
    double res;
    struct res_by_kernel gflops;
    solve(N, No, IA, x, JA, A, b, &count_steps, &res, maxiter, tol, offset_on_proc, indexes_on_send, indexes_on_recv, &gflops);

    end_time = MPI_Wtime();

    struct res_by_kernel gflops_all;
    if (rank == ROOT) {
        MPI_Reduce(&gflops.time_dot, &gflops_all.time_dot, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.time_axpby, &gflops_all.time_axpby, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.time_SpMV, &gflops_all.time_SpMV, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.time_VVbe, &gflops_all.time_VVbe, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_dot, &gflops_all.gflops_dot, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_axpby, &gflops_all.gflops_axpby, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_SpMV, &gflops_all.gflops_SpMV, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_VVbe, &gflops_all.gflops_VVbe, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(&gflops.time_dot, NULL, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.time_axpby, NULL, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.time_SpMV, NULL, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.time_VVbe, NULL, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_dot, NULL, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_axpby, NULL, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_SpMV, NULL, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
        MPI_Reduce(&gflops.gflops_VVbe, NULL, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
    }

    if (rank == ROOT) {
        printf("\nSolved\n");
        printf("Time of calculating: %.5f seconds\n", end_time - start_time);
        // printf("Exepted memory: %d Kbytes\n", exepted_memory/1024);

        printf("\nTiming by kernels\n");
        printf("Time of dot: %.5f seconds, gflops: %.5f\n", gflops_all.time_dot, gflops_all.gflops_dot/gflops_all.time_dot);
        printf("Time of axpby: %.5f seconds, gflops: %.5f\n", gflops_all.time_axpby, gflops_all.gflops_axpby/gflops_all.time_axpby);
        printf("Time of SpMV: %.5f seconds, gflops: %.5f\n", gflops_all.time_SpMV, gflops_all.gflops_SpMV/gflops_all.time_SpMV);
        printf("Time of VVbe: %.5f seconds, gflops: %.5f\n\n", gflops_all.time_VVbe, gflops_all.gflops_VVbe/gflops_all.time_VVbe);
        printf(
            "Time of solver: %.5f seconds, gflops: %.5f\n\n",
            gflops_all.time_dot+gflops_all.time_axpby+gflops_all.time_SpMV+gflops_all.time_VVbe,
            (gflops_all.gflops_dot+gflops_all.gflops_axpby+gflops_all.gflops_SpMV+gflops_all.gflops_VVbe)
            / (gflops_all.time_dot+gflops_all.time_axpby+gflops_all.time_SpMV+gflops_all.time_VVbe)
         );
    }

    MPI_Finalize();
    return 0;
}