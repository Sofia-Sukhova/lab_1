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

    output = fopen("output.csv", "w");


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
    double *u;
    u = calloc(K * M + 10, sizeof(double));

    for (m = 0; m < M; m++){
        u[m] = fi(m * x_step);
    }

    int num = M / commsize;
    MPI_Status status;
    MPI_Request recv;
    if (rank == commsize - 1){
        for (k = 0; k < K - 1; k ++){
            for (m = rank * num; m < M; m++){
                if (m == num * rank){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    u[k * M + m - 1] = tmp;
                }

                u = sceme_realization(u, k, m);

                if (m == num * rank){
                    tmp = u[k * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &recv);
                }
            }
        }  

        // match all matrix together and put to the file
        for(i = 0; i < commsize - 1; i++){
            double *uTMP = calloc(K * M + 10, sizeof(double));

            // MPI_Recv(uTMP, K * M, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(uTMP, K * num, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            printf("here i am, rank = %d\n", i);

            for (k = 0; k < K; k++){
                for(m = 0; m < M; m ++){
                    printf("%f ", uTMP[k * M + m]);
                }
                printf("\n");
            }

            for (k = 0; k < K; k++){
                for(m = 0; m < num; m++){
                    u[k * M + m + i * num] = uTMP[k * num + m];
                }
            }
            free(uTMP);
        }

        fprintf(output, "t\tx\tu\n");
        for (i = 0; i < K; i++){
            for (j = 0; j < M ; j++){
                fprintf(output, "%lf\t%lf\t%lf\n", i * t_step, x_step * j, u[i * M  + j]);
            }
        } 
    }
    if (rank == 0){
        for (k = 0; k < K; k ++){
            for (m = 0; m < num; m++){
                if ((m == num - 1) && (k >= 1)){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    u[k * M + m + 1] = tmp;
                }
                u = sceme_realization(u, k, m);
                if (m == num - 1){
                    tmp = u[k * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recv);
                }
            }
        }
        MPI_Send(u, K * M, MPI_DOUBLE, commsize - 1, 0, MPI_COMM_WORLD);
    }
    if ((rank != commsize - 1) && (rank != 0)){
        int x_end = num * (rank + 1);
        for (k = 0; k < K; k ++){
            for (m = rank * num; m < x_end; m++){
                if ((m == x_end - 1) && (k >= 1)){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    u[k * M + m + 1] = tmp;
                }
                if (m == num * rank){
                    MPI_Recv(&tmp, 1, MPI_DOUBLE, rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    u[k * M + m - 1] = tmp;
                }
                u = sceme_realization(u, k, m);
                if (m == x_end - 1){
                    tmp = u[k * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recv);
                }
                if (m == num * rank){
                    tmp = u[k * M + m];
                    MPI_Isend(&tmp, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &recv);
                }
            }
        } 
        sleep(rank);
        printf("HERE, rank = %d\n", rank);
        for (k = 0; k < K; k++){
            for (m = 0; m < M; m ++){
                printf("%f ", u[k * M + m]);
            }
            printf("\n");
        }
        // MPI_Send(u, K * M, MPI_DOUBLE, commsize - 1, 0, MPI_COMM_WORLD); 

        
        MPI_Send(u, K * num, MPI_DOUBLE, commsize - 1, 0, MPI_COMM_WORLD);     

    }
    free(u);
    fclose(output);
    MPI_Finalize();
}