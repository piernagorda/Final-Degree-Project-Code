#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

#include "blis.h"
#include "defs.h"


//void sgemm_asm_8x6(int k, double* restrict alpha, float* restrict a, float* restrict b, double* restrict beta, float* restrict c, int rs_c, int cs_c);
void sgemm_ref(int k, int mr_alg, int nr_alg, 	float* restrict alpha, float* restrict a, float* restrict b, float* restrict beta, float* restrict c, int rs_c, int cs_c);
void sgemm_armv8a_asm_8x12(dim_t m_alg, dim_t n_alg,dim_t k_alg, float* restrict alphap, float* restrictAr, float* restrict Br,float* restrict beta, float* restrictCr, inc_t rs_c, inc_t ldc);

void Pack_A(float *A, unsigned int lda, float *A_pack, unsigned int m, unsigned int k, int * threads)
{
	float *A_pack_local;

	#pragma omp parallel for num_threads(*threads) private(A_pack_local)
	for(unsigned int ic=0;ic<m;ic+=BLOCK_MR){
		
		A_pack_local=&A_pack[ic*k];
		unsigned int m_alg=fmin(BLOCK_MR,m-ic);
		for(unsigned int pc=0;pc<k;pc++){
			
			for(unsigned int ir=0;ir<m_alg;ir++){
				A_pack_local[0]=A[(ic+ir)+pc*lda];
				A_pack_local++;
			}
		}
		
	}
}

void Pack_A_trans(float *A, unsigned int lda, float *A_pack, unsigned int m, unsigned int k, int * threads)
{
        float *A_pack_local;

        //#pragma omp parallel for num_threads(*threads) private(A_pack_local)
        for(unsigned int ic=0;ic<m;ic+=BLOCK_MR){

                A_pack_local=&A_pack[ic*k];
                unsigned int m_alg=fmin(BLOCK_MR,m-ic);
                for(unsigned int pc=0;pc<k;pc++){

                        for(unsigned int ir=0;ir<m_alg;ir++){
                                 A_pack_local[0]=A[(ic+ir)*lda+pc];    // Trans
                                 A_pack_local++;
                        }
                }
        }
}


void Pack_B(float *B, unsigned int ldb, float *B_pack, unsigned int k, unsigned int n, int * threads)
{
        float *B_pack_local;

	#pragma omp parallel for num_threads(*threads) private(B_pack_local)
        for(unsigned int jc=0;jc<n;jc+=BLOCK_NR){

		B_pack_local=&B_pack[jc*k];
                unsigned int n_alg=fmin(BLOCK_NR,n-jc);
                for(unsigned int pc=0;pc<k;pc++){

                        for(unsigned int jr=0;jr<n_alg;jr++){
                                B_pack_local[0]=B[pc+jc*ldb+jr*ldb];
                                B_pack_local++;
                        }
                }

        }
}


void Pack_B_trans(float *B, unsigned int ldb, float *B_pack, unsigned int k, unsigned int n, int * threads)
{
        float *B_pack_local;

        //#pragma omp parallel for num_threads(*threads) private(B_pack_local)
        for(unsigned int jc=0;jc<n;jc+=BLOCK_NR){

                B_pack_local=&B_pack[jc*k];
                unsigned int n_alg=fmin(BLOCK_NR,n-jc);
                for(unsigned int pc=0;pc<k;pc++){

                        for(unsigned int jr=0;jr<n_alg;jr++){
                                B_pack_local[0]=B[jr+jc+pc*ldb];        // Trans
                                B_pack_local++;
                        }
                }

        }
}


void mydgemm_(char *transA_, 
		char *transB_,
		int *m_, int *n_, int *k_,
		float * alphap,
		float * A, int *lda_,
		float * B, int *ldb_,
		float * betap,
		float * C, int *ldc_,
		void * Ac_pack_v, void * Bc_pack_v,
		int * threads) {

        char transA=*transA_;
        char transB=*transB_;
		float *Ac, *Bc;
		float *Cc;
		float *Ar, *Br;
		float *Cr;
		float beta;
		unsigned int m=(unsigned int)*m_, n=(unsigned int)*n_,k=(unsigned int)*k_;
		unsigned int lda=(unsigned int)*lda_, ldb=(unsigned int)*ldb_, ldc=(unsigned int)*ldc_;
		float *Ac_pack=(float *)Ac_pack_v;
		float *Bc_pack=(float *)Bc_pack_v;

        for (unsigned int jc=0; jc<n; jc+=BLOCK_NC) {
		//printf("Entrando en bucle 1 \n");
		unsigned int n_alg=fmin(BLOCK_NC,n-jc);
            for (unsigned int pc=0; pc<k; pc+=BLOCK_KC) {
			unsigned int k_alg=fmin(BLOCK_KC,k-pc);
			if (pc >= BLOCK_KC) //Check beta
				beta=1.0;
			else beta=*betap;
			if(transB=='N'){
				Bc=&B[pc+jc*ldb];
				Pack_B(Bc, ldb, Bc_pack, k_alg, n_alg, threads);  //PACK B
			}
			else{
				Bc=&B[pc*ldb+jc];
                Pack_B_trans(Bc, ldb, Bc_pack, k_alg, n_alg, threads);  //PACK B
			}
              for (unsigned int ic=0; ic<m; ic+=BLOCK_MC) {
				unsigned int m_alg=fmin(BLOCK_MC,m-ic);
				float *Ac_pack_local=&Ac_pack[omp_get_thread_num()*BLOCK_MC*BLOCK_KC]; // Ac pack pointer per Loop 3 thread
				if(transA=='N'){
					Ac=&A[ic+pc*lda];
					Pack_A(Ac,lda,(float*)Ac_pack_local,m_alg,k_alg, threads); //PACK A
				}
				else{
					Ac=&A[ic*lda+pc];
					Pack_A_trans(Ac,lda,(float*)Ac_pack_local,m_alg,k_alg, threads); //PACK A	
				}
				Cc=&C[ic+jc*ldc];
				#pragma omp parallel num_threads(*threads) private(Ar, Br, Cr)
                                {
					#pragma omp for private(Ar, Br, Cr) 
					for(unsigned jr=0;jr<n_alg;jr+=BLOCK_NR){
						unsigned int nr_alg=fmin(BLOCK_NR,n_alg-jr);
						for(unsigned int ir=0;ir<m_alg;ir+=BLOCK_MR){
							unsigned int mr_alg=fmin(BLOCK_MR,m_alg-ir);
							Ar=&Ac_pack_local[ir*k_alg];
							Br=&Bc_pack[jr*k_alg];
							Cr=&Cc[ir+jr*ldc];
							if(mr_alg==BLOCK_MR && nr_alg==BLOCK_NR){
								sgemm_armv8a_asm_8x12(mr_alg, nr_alg, k_alg, alphap, Ar, Br, &beta, Cr, 1, ldc);
								//dgemm_asm_8x6(k_alg,alphap,Ar,Br,&beta,Cr,1,ldc);
							}
							else{//Micro-kernel cannot be applied
								sgemm_ref(k_alg,mr_alg,nr_alg,alphap,Ar,Br,&beta,Cr,1,ldc);
							}
						}
					}
				}
                        }
                }
        }

}
