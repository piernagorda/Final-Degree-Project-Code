#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>
#include "defs.h"


void dgemm_asm_8x6(int k, double* restrict alpha, double* restrict a, double* restrict b, double* restrict beta, double* restrict c, int rs_c, int cs_c);

void dgemm_ref(int k, int mr_alg, int nr_alg, 	double* restrict alpha, double* restrict a, double* restrict b, double* restrict beta, double* restrict c, int rs_c, int cs_c);


void Pack_A(double *A, unsigned int lda, double *A_pack, unsigned int m, unsigned int k, int * threads)
{
	double *A_pack_local;

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

void Pack_A_trans(double *A, unsigned int lda, double *A_pack, unsigned int m, unsigned int k, int * threads)
{
        double *A_pack_local;

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


void Pack_B(double *B, unsigned int ldb, double *B_pack, unsigned int k, unsigned int n, int * threads)
{
        double *B_pack_local;

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


void Pack_B_trans(double *B, unsigned int ldb, double *B_pack, unsigned int k, unsigned int n, int * threads)
{
        double *B_pack_local;

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


void mydgemm_(char *transA_, char *transB_,
			int *m_, int *n_, int *k_,
			double * alphap,
                        double * A, int *lda_,
                        double * B, int *ldb_,
			double * betap,
                        double * C, int *ldc_,
			void * Ac_pack_v, void * Bc_pack_v,
			int * threads) {

        char transA=*transA_;
        char transB=*transB_;
	double *Ac, *Bc;
	double *Cc;
	double *Ar, *Br;
	double *Cr;
	double beta;

	unsigned int m=(unsigned int)*m_, n=(unsigned int)*n_,k=(unsigned int)*k_;
	unsigned int lda=(unsigned int)*lda_, ldb=(unsigned int)*ldb_, ldc=(unsigned int)*ldc_;

	double *Ac_pack=(double *)Ac_pack_v;
	double *Bc_pack=(double *)Bc_pack_v;

        for (unsigned int jc=0; jc<n; jc+=BLOCK_NC) {

		unsigned int n_alg=fmin(BLOCK_NC,n-jc);
                for (unsigned int pc=0; pc<k; pc+=BLOCK_KC) {

			unsigned int k_alg=fmin(BLOCK_KC,k-pc);
			if (pc >= BLOCK_KC) //Check beta
				beta=1.0;
			else
				beta=*betap;
			if(transB=='N')
			{
				Bc=&B[pc+jc*ldb];
				Pack_B(Bc, ldb, Bc_pack, k_alg, n_alg, threads);  //PACK B
			}
			else
			{
				Bc=&B[pc*ldb+jc];
                                Pack_B_trans(Bc, ldb, Bc_pack, k_alg, n_alg, threads);  //PACK B
			}


                        for (unsigned int ic=0; ic<m; ic+=BLOCK_MC) {

				unsigned int m_alg=fmin(BLOCK_MC,m-ic);
				double *Ac_pack_local=&Ac_pack[omp_get_thread_num()*BLOCK_MC*BLOCK_KC]; // Ac pack pointer per Loop 3 thread
				
				if(transA=='N')
				{
					Ac=&A[ic+pc*lda];
					Pack_A(Ac,lda,(double*)Ac_pack_local,m_alg,k_alg, threads); //PACK A
				}
				else
				{
					Ac=&A[ic*lda+pc];
					Pack_A_trans(Ac,lda,(double*)Ac_pack_local,m_alg,k_alg, threads); //PACK A	
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

							if(mr_alg==BLOCK_MR && nr_alg==BLOCK_NR)
							{
								//dgemm_asm_8x6(k_alg,alphap,Ar,Br,&beta,Cr,1,ldc);
								dgemm_ref(k_alg,mr_alg,nr_alg,alphap,Ar,Br,&beta,Cr,1,ldc);
							}
							else{//Micro-kernel cannot be applied

								dgemm_ref(k_alg,mr_alg,nr_alg,alphap,Ar,Br,&beta,Cr,1,ldc);
							}
						}
					}
				}
                        }
                }
        }

}
