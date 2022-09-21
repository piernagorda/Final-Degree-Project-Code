all: test_gemm
	
clean:
	rm -f *.o *~ core *.x

test_gemm: matrix.o
	gcc-11 sgemm-OpenBLAS.c matrix.c -I/Users/javierp/OpenBLAS/  -L/Users/javierp/OpenBLAS -lm -lopenblas -fopenmp -fPIC -o test_gemm.x

matrix.o: matrix.c
	gcc -c matrix.c

.PHONY: all clean






