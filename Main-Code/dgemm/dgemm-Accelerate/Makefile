all: test_gemm
	
clean:
	rm -f *.o *~ core *.x

test_gemm: main.o matrix.o
	gcc-11 main.o matrix.o -framework Accelerate -flax-vector-conversions -o test_gemm.x
main.o: main.c
	gcc -c main.c
matrix.o: matrix.c
	gcc -c matrix.c

.PHONY: all clean






