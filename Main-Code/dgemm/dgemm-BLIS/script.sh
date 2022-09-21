

export BLIS_JC_NT=1 #bucle 1 (N)
export BLIS_IC_NT=1 #bucle 3 (M)
export BLIS_JR_NT=1 #bucle 4 (N)
export BLIS_IR_NT=1 #bucle 5 (M)

#export OMP_PROC_BIND=true
#export KMP_HOT_TEAMS_MODE=1
#export OMP_WAIT_POLICY=active
export OMP_NUM_THREADS=1
export OMP_NESTED=true


echo "GEMM 1"
export BLIS_JR_NT=1 #bucle 4 (N)
export OMP_NUM_THREADS=1
./test_gemm.x

echo "GEMM 2"
export BLIS_JR_NT=2 #bucle 4 (N)
export OMP_NUM_THREADS=2
./test_gemm.x

echo "GEMM 4"
export BLIS_JR_NT=4 #bucle 4 (N)
export OMP_NUM_THREADS=4
./test_gemm.x

