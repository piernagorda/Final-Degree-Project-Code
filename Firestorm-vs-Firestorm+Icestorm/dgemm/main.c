#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <Accelerate/Accelerate.h>
#include "defs.h"
#include <pthread/qos.h>
#include <pthread.h>

struct thread_args{
	void *A;
	void *B;
	void *C;
	int m;
	int n;
	int k;
	float alphap;
	float betap;
	int transa;
	int transb;
};

void random_matrix(    double *M, unsigned int height, unsigned int width);

void print_matrix(    double *M, unsigned int height, unsigned int width);

void compare_matrix(double *MG, double *MT, unsigned int height, unsigned int width);

void firestormGEMM(void * arguments){
    struct thread_args *struct_ptr = (struct thread_args*) arguments;
	int lda=struct_ptr->k, ldb=struct_ptr->k, ldc=struct_ptr->k;
	int e = pthread_set_qos_class_self_np(QOS_CLASS_USER_INTERACTIVE, 0);
	if (e) {
		fprintf(stderr, "Pthread error: %d\n", e);
		exit(1);
	}
	else printf("QoS: QOS_CLASS_USER_INTERACTIVE \n");
	printf("firestorm initiated \n");
	//sgemm_(&struct_ptr->transa, &struct_ptr->transb, &struct_ptr->m, &struct_ptr->n, &struct_ptr->k, &struct_ptr->alphap, (float*)struct_ptr->A, &lda, (float*)struct_ptr->B, &ldb, &struct_ptr->betap, (float*)struct_ptr->C, &ldc);
	cblas_dgemm(CblasColMajor, struct_ptr->transa, struct_ptr->transb, struct_ptr->m, struct_ptr->n, struct_ptr->k, struct_ptr->alphap, struct_ptr->A, lda, struct_ptr->B, ldb, struct_ptr->betap, struct_ptr->C, ldc);
    printf("firestorm finished \n");
}

void icestormGEMM(void * arguments){
	struct thread_args *struct_ptr = (struct thread_args*) arguments;
	int lda=struct_ptr->k, ldb=struct_ptr->k, ldc=struct_ptr->k;
	int e = pthread_set_qos_class_self_np(QOS_CLASS_BACKGROUND, 0);
	if (e) {
		fprintf(stderr, "Pthread error: %d\n", e);
		exit(1);
	}
	else printf("QoS: QOS_CLASS_BACKGROUND \n");
	printf("icestorm initiated \n");
	//sgemm_(&struct_ptr->transa, &struct_ptr->transb, &struct_ptr->m, &struct_ptr->n, &struct_ptr->k, &struct_ptr->alphap, (float*)struct_ptr->A, &lda, (float*)struct_ptr->B, &ldb, &struct_ptr->betap, (float*)struct_ptr->C, &ldc);
	cblas_dgemm(CblasColMajor, struct_ptr->transa, struct_ptr->transb, struct_ptr->m, struct_ptr->n, struct_ptr->k, struct_ptr->alphap, struct_ptr->A, lda, struct_ptr->B, ldb, struct_ptr->betap, struct_ptr->C, ldc);
    printf("icestorm finished \n");
}

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

            setQoS(3);
            int threads = 1;
            char str[3];
            sprintf(str, "%d", threads);
            setenv("VECLIB_MAXIMUM_THREADS", str, true);
            printf("Number of Threads : %s\n", getenv("VECLIB_MAXIMUM_THREADS"));

            void *A, *B, *C, *C_save;
            void *Ac_pack_v, *Bc_pack_v;
            assert(!posix_memalign(&Ac_pack_v, 4096, MAX_THREAD * BLOCK_MC * BLOCK_KC * sizeof(double)));
                assert(!posix_memalign(&Bc_pack_v, 4096, MAX_THREAD * BLOCK_KC * BLOCK_NC * sizeof(double)));
            printf("Results of Accelerate: \n");
            printf("\n");
            for (int dim=1000  ;dim <= 31000;dim+=5000){
                int m=dim, n=dim, k=dim;
                int lda=m, ldb=k, ldc=m;
                double GFLOPS_prueba=0.0, GFLOPS_BLIS=0.0;
                struct timeval start,end;
                double time,time_ref;
                double ops=2.0*m*n*k;

                assert(!posix_memalign(&A, 4096, m*k*sizeof(double)));
                assert(!posix_memalign(&B, 4096, k*n*sizeof(double)));
                assert(!posix_memalign(&C, 4096, m*n*sizeof(double)));
                assert(!posix_memalign(&C_save, 4096, m*n*sizeof(double)));
                random_matrix((double*)A,m,k);
                //printf("Matrix A\n");
                //print_matrix(A, m, k);
                random_matrix((double*)B,k,n);
                //printf("Matrix B\n");
                //print_matrix(B, k, n);
                random_matrix((double*)C,m,n);
                //printf("Matrix C\n");
                //print_matrix(C, m, n);
                memcpy(C_save, C, m*n*sizeof(double));
                double  alphap=2.0;
                double  betap=2.0;

                int transa = 111;
                int transb = 111;
                char transaOut = 'N';
                char transbOut = 'N';
                
                gettimeofday( &start, NULL );
                cblas_dgemm(CblasColMajor, transa, transb, m, n, k, alphap, A, lda, B, ldb, betap, C, ldc);
                gettimeofday( &end, NULL );
                time_ref=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
                GFLOPS_BLIS=ops/(time_ref*1.0e9);
                time=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;

                printf("%c %c\t%d\t%d\t%d\t%5.6f\t%3.3f\t%d\t \n",transaOut,transbOut,m,n,k,time,GFLOPS_BLIS,threads);


                //INTENTO DE PARALELIZACION

                int m_asup = (m*90)/100; //Asignamos un 92% de filas al Firestorm
                int m_ainf = dim - m_asup; //El otro 8% al Icestorm


                pthread_t th1, th2;
                void *A_inferior = &A[m_asup*sizeof(double)];
                void *C_inferior = &C_save[m_asup*sizeof(double)];

                //((float*)C_save)[m_asup+dim]=1.0;

                struct thread_args argumentosFirestorm[1];
                argumentosFirestorm[0].A = A;
                argumentosFirestorm[0].B = B;
                argumentosFirestorm[0].C = C_save;
                argumentosFirestorm[0].m=m_asup;
                argumentosFirestorm[0].n=n;
                argumentosFirestorm[0].k=k;
                argumentosFirestorm[0].alphap=alphap;
                argumentosFirestorm[0].betap=betap;
                argumentosFirestorm[0].transa=transa;
                argumentosFirestorm[0].transb=transb;

                struct thread_args argumentosIcestorm[1];
                argumentosIcestorm[0].A = A_inferior;
                argumentosIcestorm[0].B = B;
                argumentosIcestorm[0].C = C_inferior;
                argumentosIcestorm[0].m=m_ainf;
                argumentosIcestorm[0].n=n;
                argumentosIcestorm[0].k=k;
                argumentosIcestorm[0].alphap=alphap;
                argumentosIcestorm[0].betap=betap;
                argumentosIcestorm[0].transa=transa;
                argumentosIcestorm[0].transb=transb;

                gettimeofday( &start, NULL );
                pthread_create(&th1, NULL, (void*)icestormGEMM, &argumentosIcestorm);
                pthread_create(&th2, NULL, (void*)firestormGEMM, &argumentosFirestorm);
                pthread_join(th2,NULL);
                pthread_join(th1,NULL);
                gettimeofday( &end, NULL );

                time_ref=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
                GFLOPS_prueba=ops/(time_ref*1.0e9);
                printf("GFLOPS PRUEBA: %3.3f \n", GFLOPS_prueba);
                compare_matrix(C, C_save, dim, dim);
                printf("---\n");

                free(A);
                free(B);
                free(C);
                free(C_save);
             }
            free(Ac_pack_v);
            free(Bc_pack_v);
    return 0;
}
