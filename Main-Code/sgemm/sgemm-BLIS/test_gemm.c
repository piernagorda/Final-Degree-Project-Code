#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <pthread/qos.h>
#include <omp.h>
#include "blis.h"
#include "defs.h"

void mydgemm_(char* transA, char* transB,
			int* m, int* n, int* k,
			float * alphap,
			float * A, int* lda,
			float * B, int* ldb,
			float * betap,
			float * C, int* ldc,
			void * Ac_pack_v, void * Bc_pack_v,
			int * threads);

void compare_matrix(	float *MG, float *MT, unsigned int height, unsigned int width);

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
		void *A, *B, *C, *C_save;
		void *Ac_pack_v, *Bc_pack_v;
		assert(!posix_memalign(&Ac_pack_v, 4096, MAX_THREAD * BLOCK_MC * BLOCK_KC * sizeof(float)));
		assert(!posix_memalign(&Bc_pack_v, 4096, MAX_THREAD * BLOCK_KC * BLOCK_NC * sizeof(float)));

		for (int dim=1000;dim <= 31000;dim+=5000){
			int m=dim, n=dim, k=dim;
			int lda=m, ldb=k, ldc=m;
			double GFLOPS=0.0, GFLOPS_BLIS=0.0;
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
			float  alphap=2.0;
			float  betap=2.0;

			f77_char transa = 'N';
			f77_char transb = 'N';
			//Obtiene la hora de inicio
			gettimeofday( &start, NULL );
			//Procede a la multiplicacion
			sgemm_( &transa, &transb, &m, &n, &k, &alphap, (float*)A, &lda, (float*)B, &ldb, &betap, (float*)C, &ldc );
			//Obtiene hora de fin
			gettimeofday( &end, NULL );
			printf("Hecho el sgemm de BLIS \n");
			//Calcula los GFLOPS
			time_ref=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
			GFLOPS_BLIS=ops/(time_ref*1.0e9);
			//printf("%d\t%d\t%d\t%5.6f\t%3.3f\n",m,n,k,time_ref,GFLOPS_BLIS);
			//printf(" \n Matrix Out de REFERENCIA \n");
			//print_matrix(C, m, n);
			
			int threads=4;
			char *s = getenv("BLIS_JR_NT");
			if( s != NULL ) {
					sscanf(s,"%d",&threads);
			}
			
			gettimeofday( &start, NULL );
			//Hace la misma multiplicacion con la version mejorada de BLIS
			mydgemm_(&transa, &transb, &m, &n, &k, &alphap, (float*)A, &lda, (float*)B, &ldb, &betap, (float*)C_save, &ldc, Ac_pack_v, Bc_pack_v,&threads);
			gettimeofday( &end, NULL );
			printf("Hecho el sgemm de BLIS adaptado \n");
			time=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
			GFLOPS=ops/(time*1.0e9);
			//N N	1000	1000	1000	0.331291	37.281	6.037	0.16	Compare OK.
			printf("%c %c\t%d\t%d\t%d\t%5.6f\t%3.3f\t%3.3f\t%d\t",transa,transb,m,n,k,time,GFLOPS_BLIS,GFLOPS,threads);
			//printf(" \n Matrix Out de SGEMM \n");
			//print_matrix(C_save, m, n);
			//compare_matrix(C_save, C, m, n);	
			printf("---\n");
			free(A);	
			free(B);	
			free(C);	
			free(C_save);	
		}
		free(Ac_pack_v);
		free(Bc_pack_v);
	}
	return 0;
}
