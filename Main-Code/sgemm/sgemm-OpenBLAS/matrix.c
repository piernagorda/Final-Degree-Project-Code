#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>

void compare_matrix(double *MG, double *MT, unsigned int height, unsigned int width) {

	int e=0;
	for(unsigned int x=0 ; x < width ; x++){
                for(unsigned y=0 ; y < height ; y++){
                        double golden = MG[x*height+y];
			double test = MT[x*height+y];

			double delta = (golden-test);
			if (fabs(delta) > (fabs(golden) / 3000000)) {
				printf("Compare error: x=%d, y=%d, golden=%g, test=%g, delta=%g\n",x,y,golden,test,delta);
				e++;
				if (e>20) {
					printf("Too many errors, not printing any more.\n");
					return;
				}
			}

                }
        }
	if(e==0)
		printf("Compare OK.\n");
}


void random_matrix(float *M, unsigned int height, unsigned int width) {

	for(unsigned int x=0 ; x < width ; x++){
		for(unsigned y=0 ; y < height ; y++){
			M[x*height+y]= ((float) (rand() % 10000)) / 10000.0f;
		}
	}
}

void random_matrix_tri_lower(float *M, unsigned int height, unsigned int width) {

	for(unsigned int x=0 ; x < width ; x++){
		for(unsigned y=x ; y < height ; y++){
			if (y==x )	
				M[x*height+y]= 1.0 ;
			else
				M[x*height+y]= ((float) (rand() % 10000)) / 10000.0f;
		}
	}
}



void print_matrix(float *M, unsigned int height, unsigned int width) {
    for (unsigned int y=0 ; y<height ; y++) {
        for (unsigned int x=0 ; x<width ; x++) {
            printf("%f\t",M[x*height+y]);
        }
        printf("\n");
    }
}

