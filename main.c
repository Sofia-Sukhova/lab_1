#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
// #include "func.c"

extern double func  (double t, double x);
extern double fi    (double x);
extern double ksi   (double t);
extern double *sceme_realization (double *u, int k, int m);

extern double t_max;
extern double x_max;
extern double t_step;
extern double x_step;


int main( int argc, char **argv ){

    FILE *output;
    int rank, commsize, stat;
    unsigned long long i, j, k, m;
    int K = (int)(t_max/t_step);
    int M = (int)(x_max/x_step);
    double tmp = 0;
    double time_start, time_finish;

    output = fopen("output.csv", "w");


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    double *u;
    u = calloc(K * M, sizeof(double));

    for (m = 0; m < M; m++){
        u[m] = fi(m * x_step);
    }

    int num = M / commsize;
    MPI_Status status;
    MPI_Request recv;

    time_start = MPI_Wtime();

    if (rank == commsize - 1){
        for (k = 0; k < K - 1; k ++){
            for (m = rank * num; m < M; m++){
                if (commsize > 1){
                    if (m == num * rank){
                        MPI_Recv(&tmp, 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, &status);
                        u[k * M + m - 1] = tmp;
                    }

                    u = sceme_realization(u, k, m);

                    if ((m == num * rank)&&(k != K - 1)){
                        tmp = u[(k + 1) * M + m];
                        MPI_Isend(&tmp, 1, MPI_DOUBLE, rank - 1, k + 1, MPI_COMM_WORLD, &recv);
                    }
                } else{
                    u = sceme_realization(u, k, m);
                }
            }
        }
        
        time_finish = MPI_Wtime();

        // match all matrix together and put to the file
        for(i = 0; i < commsize - 1; i++){
            double *uTMP = calloc(K * M, sizeof(double));

            MPI_Recv(uTMP, K * M, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);

            for (k = 0; k < K; k++){
                for(m = i * num; m < (i + 1)*num; m++){
                    u[k * M + m] = uTMP[k * M + m];
                }
            }

            free(uTMP);
        }

        printf("start = %f, finish= %f, total = %f\n", time_start, time_finish, time_finish - time_start);

        fprintf(output, "t\tx\tu\n");
        fprintf(output, "%d\t%d\t%d\n", K, M, 0);
        for (i = 0; i < K; i++){
            for (j = 0; j < M ; j++){
                fprintf(output, "%lf\t%lf\t%lf\n", i * t_step, x_step * j, u[i * M  + j]);
            }
        } 

    }

    if ((rank == 0) && (commsize > 1)){
        for (k = 0; k < K; k ++){
            for (m = 0; m < num; m++){
                if ((m == num - 1) && (k >= 1)){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, &status);
                    u[k * M + m + 1] = tmp;
                }
                u = sceme_realization(u, k, m);
                if (m == num - 1){
                    tmp = u[k * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, &recv);
                }
            }
        }
        MPI_Send(u, K * M, MPI_DOUBLE, commsize - 1, 1, MPI_COMM_WORLD);
    }
    if ((rank != commsize - 1) && (rank != 0)){
        int x_end = num * (rank + 1);
        for (k = 0; k < K; k ++){
            for (m = rank * num; m < x_end; m++){
                if ((m == x_end - 1) && (k >= 1)){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, &status);
                    u[k * M + m + 1] = tmp;
                }
                if (m == num * rank){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank - 1, k, MPI_COMM_WORLD, &status);
                    u[k * M + m - 1] = tmp;
                }

                u = sceme_realization(u, k, m);

                if (m == x_end - 1){
                    tmp = u[k * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank + 1, k, MPI_COMM_WORLD, &recv);
                }
                if ((m == num * rank) && (k != K - 1)){
                    tmp = u[(k + 1) * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank - 1, k + 1, MPI_COMM_WORLD, &recv);
                }
            }
        } 

        MPI_Send(u, K * M, MPI_DOUBLE, commsize - 1, 1, MPI_COMM_WORLD);     

    }
    free(u);
    fclose(output);
    MPI_Finalize();
}
