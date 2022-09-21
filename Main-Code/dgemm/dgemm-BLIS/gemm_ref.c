

#include <stdio.h>
#include <stdint.h>



void dgemm_ref( int k, int mr_alg, int nr_alg, double* restrict alpha, double* restrict a, double* restrict b, double* restrict beta, double* restrict c, int rs_c, int cs_c)
{
        double ab[mr_alg*nr_alg];
        double bj, ai;
        unsigned int l,i,j;
        for ( i = 0; i < mr_alg * nr_alg ; ++i ){ //set 0s

                *(ab+i)=0;
        }

        for ( l = 0; l < k; ++l ){
                double *abij=ab;
                for ( j = 0; j < nr_alg; ++j ){
                        bj = *(b + j);
                        for ( i = 0; i < mr_alg; ++i ){
                                ai = *(a + i);
                                abij[0] += ai * bj;
                                abij++;
                        }
                }
                a+=mr_alg;
                b+=nr_alg;
        }

        for ( i = 0; i < mr_alg * nr_alg; ++i ){ //Scale by alpha

                ab[i]*=*alpha;
        }

        double *Cr_l;
        double *ab_l=ab;

        if (beta ==0){ // no C loaded from memory

                for ( j = 0; j < nr_alg; ++j ){
                        Cr_l=&c[j*cs_c];
                        for ( i = 0; i < mr_alg; ++i ){
                                *Cr_l=*ab_l;
                                ab_l++;
                                Cr_l+=rs_c;
                        }
                }
        }
        else{ // C loaded from memory and scaled by beta

                for ( j = 0; j < nr_alg; ++j ){
                        Cr_l=&c[j*cs_c];
                        for ( i = 0; i < mr_alg; ++i ){
                                *Cr_l*=*beta;
                                *Cr_l+=*ab_l;
                                ab_l++;
                                Cr_l+=rs_c;
                        }
                }
        }
}



