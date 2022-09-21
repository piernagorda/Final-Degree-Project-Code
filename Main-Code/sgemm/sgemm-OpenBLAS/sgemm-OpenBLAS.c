#include <stdio.h>
#include "cblas.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <pthread/qos.h>
#include <omp.h>
#include "defs.h"

void random_matrix(	float *M, unsigned int height, unsigned int width);

void print_matrix(	float *M, unsigned int height, unsigned int width);

void setQoS(int n){
    if (n==0){
        int e = pthread_set_qos_class_self_np(QOS_CLASS_BACKGROUND, 0);
            if (e) {
                fprintf(stderr, "Pthread error: %d\n", e);
                exit(1);
            }
            else printf("QoS: QOS_CLASS_BACKGROUND \n");
    }
    if(n==1){
        int e = pthread_set_qos_class_self_np(QOS_CLASS_UTILITY, 0);
        if (e) {
            fprintf(stderr, "Pthread error: %d\n", e);
            exit(1);
        }
        else printf("QoS: QOS_CLASS_UTILITY \n");
    }
    if(n==2){
        int e = pthread_set_qos_class_self_np(QOS_CLASS_USER_INITIATED, 0);
        if (e) {
            fprintf(stderr, "Pthread error: %d\n", e);
            exit(1);
        }
        else printf("QoS: QOS_CLASS_USER_INITIATED \n");
    }
    if (n==3){
        int e = pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
        if (e) {
            fprintf(stderr, "Pthread error: %d\n", e);
            exit(1);
        }
        else printf("QoS: QOS_CLASS_USER_INTERACTIVE \n");
    }
}


int main(int argc, char **argv) {
	for (int qos = 0; qos < 4;++qos){
        printf("\n");
        setQoS(qos);
		for (int threads=1; threads <=4; threads=threads*2){
			
			void *A, *B, *C, *C_save;
			void *Ac_pack_v, *Bc_pack_v;
			assert(!posix_memalign(&Ac_pack_v, 4096, MAX_THREAD * BLOCK_MC * BLOCK_KC * sizeof(float)));
			assert(!posix_memalign(&Bc_pack_v, 4096, MAX_THREAD * BLOCK_KC * BLOCK_NC * sizeof(float)));
			
			for (int dim=1000;dim <= 31000;dim+=5000){
				int m=dim, n=dim, k=dim;
				int lda=m, ldb=k, ldc=m;
				openblas_set_num_threads(threads);
				double GFLOPS=0.0, GFLOPS_OpenBlas=0.0;
				struct timeval start,end;
				double time,time_ref;
				double ops=2.0*m*n*k;
				
				assert(!posix_memalign(&A, 4096, m*k*sizeof(float)));
				assert(!posix_memalign(&B, 4096, k*n*sizeof(float)));
				assert(!posix_memalign(&C, 4096, m*n*sizeof(float)));
				assert(!posix_memalign(&C_save, 4096, m*n*sizeof(float)));
				//Crea una matriz aleatoria A
				random_matrix((float*)A,m,k);
				//printf("Matrix A\n");
				//print_matrix(A, m, k);
				//Crea una matriz aleatoria B
				random_matrix((float*)B,k,n);
				//printf("Matrix B\n");
				//print_matrix(B, k, n);
				//Crea una matriz aleatoria C
				random_matrix((float*)C,m,n);
				//printf("Matrix C\n");
				//print_matrix(C, m, n);
				//Reserva en memoria el tamaño de la matriz C final
				memcpy(C_save, C, m*n*sizeof(float));
				double  alphap=2.0;
				double  betap=2.0;
				char transa = 'N';
				char transb = 'N';

				//-------------------------------GENERICO---------------------------
				//gettimeofday( &start, NULL );
				//mydgemm_(&transa, &transb, &m, &n, &k, &alphap, (double*)A, &lda, (double*)B, &ldb, &betap, (double*)C, &ldc, Ac_pack_v, Bc_pack_v,&threads);
				//gettimeofday( &end, NULL );
				//time_ref=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
				//GFLOPS=ops/(time_ref*1.0e9);
				//printf(" \n Matrix Out de mydgemm \n");
				//print_matrix(C, m, n);
				
				//-------------------------------OpenBLAS---------------------------
				gettimeofday( &start, NULL );
				/*
				typedef enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102} CBLAS_ORDER;
				typedef enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113, CblasConjNoTrans=114} CBLAS_TRANSPOSE;  
				*/
				//cblas_dgemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N, OPENBLAS_CONST blasint K,
				//		 OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc);
				//cblas_dgemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA, CBLAS_TRANSPOSE TransB, m, n, k, double alpha, A, lda, B, ldb, beta, c, ldc );
				cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, alphap, (float*)A, lda, (float*)B, ldb, betap, (float*)C_save, ldc );
				//dgemm_( &transa, &transb, &m, &n, &k, &alphap, (double*)A, &lda, (double*)B, &ldb, &betap, (double*)C, &ldc );
				gettimeofday( &end, NULL );
				time=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
				GFLOPS_OpenBlas=ops/(time*1.0e9);
				//printf(" \n Matrix Out de REFERENCIA \n");
				//print_matrix(C, m, n);

				//N N	1000	1000	1000	0.331291	37.281	6.037	0.16	Compare OK.
				
				printf("%c %c\t%d\t%d\t%d\t%5.6f\t%3.3f\t %d\t\n",transa,transb,m,n,k,time,GFLOPS_OpenBlas, threads);
				
				//printf(" \n Matrix Out de OpenBLAS \n");
				//print_matrix(C_save, m, n);
				//compare_matrix(C_save, C, m, n);	

				free(A);	
				free(B);	
				free(C);	
				free(C_save);	
				}
			free(Ac_pack_v);
			free(Bc_pack_v);
			int myInt;
			scanf("%d", &myInt);
		}
	}
	return 0;
}
