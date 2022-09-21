

#include <stdio.h>
#include <stdint.h>



void sgemm_ref( int k, int mr_alg, int nr_alg, float* restrict alpha, float* restrict a, float* restrict b, float* restrict beta, float* restrict c, int rs_c, int cs_c)
{
        float ab[mr_alg*nr_alg];
        float bj, ai;
        unsigned int l,i,j;
        for ( i = 0; i < mr_alg * nr_alg ; ++i ){ //set 0s

                *(ab+i)=0;
        }

        for ( l = 0; l < k; ++l ){
                float *abij=ab;
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

        float *Cr_l;
        float *ab_l=ab;

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



