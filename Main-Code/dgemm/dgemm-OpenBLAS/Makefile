all: test_gemm
	
clean:
	rm -f *.o *~ core *.x

test_gemm: main.o matrix.o
	gcc-11 dgemm-OpenBLAS.c matrix.c gemm.c gemm_ref.c -I//Users/javierp/OpenBLAS/  -L/Users/javierp/OpenBLAS -lm -lopenblas -fopenmp -fPIC -o test_gemm.x
main.o: main.c
	gcc -c main.c
matrix.o: matrix.c
	gcc -c matrix.c

.PHONY: all clean






