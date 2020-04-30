# CIBENCH

Welcome to the world of CIBENCH!

Enjoy the collection of 12 benchmarks and their optimization methods, they are:

1. Bzip2-1.0.6
2. Backprop, 
3. Hotspot3D,
4. LavaMD,
5. Srad_v2 from Rodinia Benchmark Suite
6. FFT from GNU scientific library
7. msb(msgrate) from NERSC8 Trinity Benchmarks
8. USQCD Chroma
9. Hmmer, 
10. H264ref,
11. Povray from SPEC CPU2006 Benchmark Suite
12. Hoard 

## Installation:

Above benchmarks (besides FFT and Chroma) can be compiled individually in their folders, or you can compile them together by following commands:

To compile the original benchmarks:
```sh
$ make org
```

To compile the optimized benchmarks:
```sh
$ make opt
```

You can easily adjust your compiler preferences by modifying $(COMPILER_CC) and $(COMPILER_CXX) in Makefile. The produced binaries enable the highest optimization levels, including -O3, -Ofast, -march=native, -mtune=native, profile guided optimization (PGO), and link timeoptimization (LTO), the optimization levels can be downgraded by modifying $(CFLAGS). Since enabling PGO optimization, the "make" process will start with an instrumented compilation, follow by a profiled execution, the information from the profiled execution will be fed back to the compiler.

To install USQCD Chroma and GSL, please following the guidance inside corresponding folders.

Due to the copyright of SPEC Benchmark Suite, we shall not release the code of Hmmer, H264ref and Povray here, please refer to our paper "What Every Scientific Programmer Should Know AboutCompiler Optimizations?" for the detailed optimization methods.

Here we list the inputs we used for your references:

For hmmer:
```sh
$hmmer: ./hmmer nph3.hmm swiss41
```
For h264ref:
```sh
$h264ref: ./h264ref -d foreman_ref_encoder_baseline.cfg
```
For povray:
```sh
$povray: ./povray SPEC-benchmark-ref.ini
```

If you have insights or questions towards CIBENCH please send email to jtan02@email.wm.edu, thanks!

