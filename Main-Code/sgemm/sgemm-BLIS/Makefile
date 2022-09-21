CFLAGS  = -O3 -march=native  -Wall -fPIC  -Wno-unused-function -std=c99 -m64 -fopenmp -D_POSIX_C_SOURCE=200112L -I/Users/javierp/local/include/blis/ -I/Users/javierp/blis/kernels/armv8a/3/ -I/Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/usr/include/pthread

LDFLAGS = -L/Users/javierp/local/lib -lblis -lm

CC	= gcc-11

all: test_gemm

clean:
	rm -f *.o *~ core *.x

test_gemm: gemm.o test_gemm.o matrix.o gemm_armv8a_asm_d6x8.o gemm_ref.o
	${CC} -fopenmp -o test_gemm.x $^ ${LDFLAGS}

.PHONY: all clean






