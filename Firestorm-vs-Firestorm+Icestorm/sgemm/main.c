#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <pthread/qos.h>
#include <Accelerate/Accelerate.h>
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

void random_matrix(	float *M, unsigned int height, unsigned int width);

void random_matrixC(float *M, unsigned int height, unsigned int width);

void print_matrix(	float *M, unsigned int height, unsigned int width);

void compare_matrix(float *MG, float *MT, unsigned int height, unsigned int width);

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
	cblas_sgemm(CblasColMajor, struct_ptr->transa, struct_ptr->transb, struct_ptr->m, struct_ptr->n, struct_ptr->k, struct_ptr->alphap, struct_ptr->A, lda, struct_ptr->B, ldb, struct_ptr->betap, struct_ptr->C, ldc);
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
	cblas_sgemm(CblasColMajor, struct_ptr->transa, struct_ptr->transb, struct_ptr->m, struct_ptr->n, struct_ptr->k, struct_ptr->alphap, struct_ptr->A, lda, struct_ptr->B, ldb, struct_ptr->betap, struct_ptr->C, ldc);
    printf("icestorm finished \n");
}

int main(int argc, char **argv) {
		setQoS(3);
		void *A, *B, *C, *C_save;
		for (int dim=1000;dim <=31000;dim+=5000){

			int m=dim, n=dim, k=dim;
			int lda=m, ldb=k, ldc=m;
			double GFLOPS_BLIS=0.0, GFLOPS_prueba=0.0;
			struct timeval start,end;
			double time_ref;
			double ops=2.0*m*n*k;

			assert(!posix_memalign(&A, 4096, m*k*sizeof(float)));
			assert(!posix_memalign(&B, 4096, k*n*sizeof(float)));
			assert(!posix_memalign(&C, 4096, m*n*sizeof(float)));
			assert(!posix_memalign(&C_save, 4096, m*n*sizeof(float)));

			random_matrix((float*)A,m,k);
			//printf("Matrix A\n");
			//print_matrix(A, m, k);
			random_matrix((float*)B,k,n);
			//printf("Matrix B\n");
			//print_matrix(B, k, n);
			random_matrix((float*)C,m,n);
			//printf("Matrix C\n");
			//print_matrix(C, m, n);
			memcpy(C_save, C, m*n*sizeof(float));
			//print_matrix(C_save, m, n);
			float  alphap=2.0;
			float  betap=2.0;

			int transa = 111;
            int transb = 111;

			gettimeofday( &start, NULL );
			cblas_sgemm(CblasColMajor, transa, transb, m, n, k, alphap, A, lda, B, ldb, betap, C, ldc);
			gettimeofday( &end, NULL );
			time_ref=(double)(((end.tv_sec-start.tv_sec) * 1000000 + (end.tv_usec-start.tv_usec)))/1000000;
			GFLOPS_BLIS=ops/(time_ref*1.0e9);
			
			printf("Accelerate %c %c\t%d\t%d\t%d\t%5.6f\t%3.3f\t",transa,transb,m,n,k,time_ref,GFLOPS_BLIS);
			//printf("Matriz C ref: \n");
			//print_matrix(C, dim, dim);
			//compare_matrix(C_save, C, m, n);	
			printf("\n");




			//INTENTO DE PARALELIZACION

			int m_asup = (m*90)/100; //Asignamos un 92% de filas al Firestorm
			int m_ainf = dim - m_asup; //El otro 8% al Icestorm


			pthread_t th1, th2;
			void *A_inferior = &A[m_asup*sizeof(float)];
			void *C_inferior = &C_save[m_asup*sizeof(float)];

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
			//printf("Matriz C paralel \n");
			//print_matrix(C_save, dim, dim);

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
	
	return 0;
}
