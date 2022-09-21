#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <Accelerate/Accelerate.h>
#include "defs.h"
#include <pthread/qos.h>

void random_matrix(    double *M, unsigned int height, unsigned int width);

void print_matrix(    double *M, unsigned int height, unsigned int width);

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
        for (int threads = 1; threads <=4; threads=threads*2){
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
            for (int dim=1000;dim <= 31000;dim+=5000){
                int m=dim, n=dim, k=dim;
                int lda=m, ldb=k, ldc=m;
                double GFLOPS=0.0, GFLOPS_BLIS=0.0;
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
                //gemm_( &transa, &transb, &m, &n, &k, &alphap, (double*)A, &lda, (double*)B, &ldb, &betap, (double*)C, &ldc );
                cblas_dgemm(CblasRowMajor, transa, transb, m, n, k, alphap, A, lda, B, ldb, betap, C, ldc);
                gettimeofday( &end, NULL );
                time_ref=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
                GFLOPS_BLIS=ops/(time_ref*1.0e9);
                //printf("%d\t%d\t%d\t%5.6f\t%3.3f\n",m,n,k,time_ref,GFLOPS_BLIS);
                //printf("Matrix Out\n");
                //print_matrix(C, m, n);

                time=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
                printf("%c %c\t%d\t%d\t%d\t%5.6f\t%3.3f\t%d\t \n",transaOut,transbOut,m,n,k,time,GFLOPS_BLIS,threads);
                //printf("Matrix Out\n");
                //print_matrix(C_save, m, n);
                //compare_matrix(C_save, C, m, n);
                free(A);
                free(B);
                free(C);
                free(C_save);
             }
            free(Ac_pack_v);
            free(Bc_pack_v);
        }
    }
    return 0;
}
