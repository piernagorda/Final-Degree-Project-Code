CFLAGS  = -O3 -march=native  -Wall -fPIC  -Wno-unused-function -std=c99 -m64 -fopenmp -D_POSIX_C_SOURCE=200112L -I/Users/javierp/local/include/blis/ -I/Users/javierp/blis/kernels/armv8a/3/

LDFLAGS = -L/Users/javierp/local/lib -lblis -lm

CC	= gcc-11

all: test_gemm

clean:
	rm -f *.o *~ core *.x

test_gemm: gemm.o test_gemm.o matrix.o gemm_armv8a_asm_d6x8.o gemm_ref.o
	${CC} -fopenmp -o test_gemm.x $^ ${LDFLAGS}


bli_gemm_aplmx_mac_s_sh32x32.o: bli_gemm_aplmx_mac_s_sh32x32.c
	gcc-11 -c bli_gemm_aplmx_mac_s_sh32x32.c  -I/Users/javierp/blis_apple-amx-dev/include/aaplmx/ -I/Users/javierp/blis_apple-amx-dev/kernels/aaplmx/3  /Users/javierp/blis_apple-amx-dev/lib/aaplmx/libblis.a -lm

.PHONY: all clean






