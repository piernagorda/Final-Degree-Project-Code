# Final Degree Project Code

In this folder you can find the following code:

- FILE: Accelerate-vs-GPU.swift
- FOLDER: Firestorm-vs-Firestorm+Icestorm
- FOLDER: Regular-Code
## Accelerate-vs-GPU.swift

This file contains a quick test for macOS systems (not M1-family-only devices) to test the performance of the Accelerate framework against the integrated GPU. It's only written for single precision (FP32) and tests the performance of the BLAS routine gemm (this is, the operation performed is sgemm). 

The implementation with Accelerate is quite straight forward, needing just three random matrixes (this is actually the slowest process) to be created and passing them as arguments to the sgemm routine. 

When it comes to implement the GPU code, the MPSMatrixMultiplication library has to be used (it's a sublibrary of Metal) and this is actually a tedious task since there's not much information available regarding this specific topic. I found [this video](https://www.youtube.com/watch?v=9pOD6wtpNWM) very useful to understand the process that has to be followed in order to send commands to the GPU.

This test can be performed in any macOS system and can even be used for iOS systems that have been jailbroken (not something common these days). It can be improved linking a C file that generates the random matrixes quicker. However, at the time of creation I didn't have enough time to investigate how to do it. The current implementation is the most efficient way to do it in Swift (which is a very slow language).

## Firestorm-vs-Firestorm+Icestorm

This folder contains a quick test for macOS systems. M1 macs (and some iPhone SoCs) have a CPU divided into two sections: the high performance cores (Firestorm cores, used for demanding workloads such as running a game, a simulator, etc) and the high efficiency cores (Icestorm cores, used for light workloads such as mail, web browsing, etc). 

The goal of this test is to see how much performance we get using one Firestorm core and then compare it against the performance that we get with a Firestorm core combined with an extra Icestorm core. 

Since the Icestorm core is highly efficient, we get an extra 7-8% increase of performance while only consuming 2-3% more of energy (this has been measured with the powermetrics command). The test can be performed on Intel Macs but doesn't really make sense since Intel CPU's are usually not divided into high performance vs high efficiency cores.

Inside the folder we can find both the single (sgemm) and double (dgemm) precision tests.

## Main-Code

This folder contains most of the code of the project. It contains tests with 1, 2 and 4 threads of OpenBLAS, BLIS and BLIS for M1 Macs (called BLIS AMX from the AMX Coprocessor). Accelerate usually runs single threaded and the number of threads cannot be modified (this is done by the GCD of macOS). When the workload is small-medium, Accelerate runs single threaded. However, when the workload increases, macOS uses more threads (a number that cannot be determined before running the tests, since it varies for each execution).

For each library we have tests in single and double precision and most of them include a Makefile. The installation of BLIS is very straight forward. However, OpenBLAS and BLIS AMX require a bit more of extra work (multithreading has to be enabled while configuring the installation and attention has to be paid!).

Sources: 

-  Accelerate: comes with macOS
-  OpenBLAS: https://www.openblas.net
-  BLIS: https://github.com/flame/blis
-  BLIS for M1: https://github.com/xrq-phys/blis_apple

