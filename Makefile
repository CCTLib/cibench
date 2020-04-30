# Choose your compiler
COMPILER_CC = gcc
COMPILER_CXX = g++

# Compilation flags for gcc to enable PGO and LTO optimizations
gcc_pl_gen = -flto -fprofile-generate
gcc_pl_use = -flto -fprofile-use -fprofile-correction
GCC_FLAGS1 = $(gcc_pl_gen)
GCC_FLAGS2 = $(gcc_pl_use)

# Compilation flags for icc to enable PGO and LTO optimizations
icc_pl_gen = -flto -prof-gen
icc_pl_use = -flto -prof-use
ICC_FLAGS1 = $(icc_pl_gen)
ICC_FLAGS2 = $(icc_pl_use)

# Compilation flags for llvm to enable PGO and LTO optimizetions
llvm_pl_gen = -flto=thin -fuse-ld=gold -fprofile-generate
llvm_pl_use = -flto=thin -fuse-ld=gold -fprofile-use=-default.profdata
LLVM_FLAGS1 = $(llvm_pl_gen)
LLVM_FLAGS2 = $(llvm_pl_use)

# Change CFLAGS1, CFLAGS2 and HFLAGS1, HFLAGS2 if you desire to use icc or llvm 
uni = -g -O3 -Ofast -march=native -mtune=native
omp = -fopenmp
CFLAGS1 = "$(uni) $(omp) $(GCC_FLAGS1)" #$(ICC_FLAGS1)/$(LLVM_FLAGS1)
CFLAGS2 = "$(uni) $(omp) $(GCC_FLAGS2)" #$(ICC_FLAGS2)/$(LLVM_FLAGS2)
HFLAGS1 = "$(uni) $(GCC_FLAGS1)"
HFLAGS2 = "$(uni) $(GCC_FLAGS2)"

# To run msb benchmarks, please have a valid MPI compiler installed and set following paths
MPICC = 
MPICC_FLAG = -cc=gcc
MPICC_INCLUDE = 
MPICC_LIBRARY = 
MPI = "$(MPICC) $(MPICC_FLAG)"
ifndef MPICC
$(error MPICC is not set)
endif

# make all original benchmarks
org: backprop_org lavaMD_org hotspot3D_org srad_org bzip_org msb_org hoard_org

# make all optimized benchmarks
opt: backprop_opt lavaMD_opt hotspot3D_opt srad_opt bzip_opt msb_op1 msb_op2 hoard_opt

# Rodinia backprop
backprop_org:
	cd ./rodinia_3.1/openmp/backprop/org && $(MAKE) CC=$(COMPILER_CC) CC_FLAGS_1=$(CFLAGS1) CC_FLAGS_2=$(CFLAGS2)
backprop_opt:
	cd ./rodinia_3.1/openmp/backprop/opt && $(MAKE) CC=$(COMPILER_CC) CC_FLAGS_1=$(CFLAGS1) CC_FLAGS_2=$(CFLAGS2)

# Rodinia lavaMD
lavaMD_org:
	cd ./rodinia_3.1/openmp/lavaMD/org && $(MAKE) CC=$(COMPILER_CC) C_C_1=$(CFLAGS1) C_C_2=$(CFLAGS2)
lavaMD_opt:
	cd ./rodinia_3.1/openmp/lavaMD/opt && $(MAKE) CC=$(COMPILER_CC) C_C_1=$(CFLAGS1) C_C_2=$(CFLAGS2)

# Rodinia hotspot3D
hotspot3D_org:
	cd ./rodinia_3.1/openmp/hotspot3D/org && $(MAKE) CC=$(COMPILER_CC) CCFLAGS1=$(CFLAGS1) CCFLAGS2=$(CFLAGS2)
hotspot3D_opt:
	cd ./rodinia_3.1/openmp/hotspot3D/opt && $(MAKE) CC=$(COMPILER_CC) CCFLAGS1=$(CFLAGS1) CCFLAGS2=$(CFLAGS2)

# Rodinia srad_v2
srad_org:
	cd ./rodinia_3.1/openmp/srad_v2/org && $(MAKE) CC=$(COMPILER_CXX) CCFLAGS1=$(CFLAGS1) CCFLAGS2=$(CFLAGS2)
srad_opt:
	cd ./rodinia_3.1/openmp/srad_v2/opt && $(MAKE) CC=$(COMPILER_CXX) CCFLAGS1=$(CFLAGS1) CCFLAGS2=$(CFLAGS2)

# bzip2
bzip_org:
	cd ./bzip2/org && $(MAKE) CC=$(COMPILER_CC) CFLAGS1=$(HFLAGS1) CFLAGS2=$(HFLAGS2)
bzip_opt:
	cd ./bzip2/opt && $(MAKE) CC=$(COMPILER_CC) CFLAGS1=$(HFLAGS1) CFLAGS2=$(HFLAGS2)

# NERSC8
# to change the wrapper compiler, i.e. from gcc to icc, please change -cc option in each Makefile
msb_org:
	cd ./NERSC/smb/org/msgrate && $(MAKE) CC=$(MPI) MPICC_INCLUDE=$(MPICC_INCLUDE) MPICC_LIBRARY=$(MPICC_LIBRARY)
msb_op1:
	cd ./NERSC/smb/opt-msb1/msgrate && $(MAKE) CC=$(MPI) MPICC_INCLUDE=$(MPICC_INCLUDE) MPICC_LIBRARY=$(MPICC_LIBRARY)
msb_opt2:
	cd ./NERSC/smb/opt-msb2/msgrate && $(MAKE) CC=$(MPI) MPICC_INCLUDE=$(MPICC_INCLUDE) MPICC_LIBRARY=$(MPICC_LIBRARY)

# Hoard
# For hoard benchmark, please adjust the compile flags here:
hoard_org:
	cd ./hoard/original/benchmarks/larson && $(MAKE) CXX=$(COMPILER_CXX) PGO_OPT1=$(HFLAGS1) PGO_OPT2=$(HFLAGS2)
hoard_opt:
	#cd ./hoard/optimized/benchmarks/larson && $(MAKE) CXX=$(COMPILER_CXX) PGO_OPT1="-flto -fprofile-generate -g -O3 -Ofast -march=native -mtune=native" PGO_OPT2="-g -O3 -Ofast -march=native -mtune=native -flto -fprofile-use -fprofile-correction"
	cd ./hoard/optimized/benchmarks/larson && $(MAKE) CXX=$(COMPILER_CXX) PGO_OPT1=$(HFLAGS1) PGO_OPT2=$(HFLAGS2)

clean:
	cd ./rodinia_3.1/openmp/backprop/org && rm -rf *.o *~ *.gcda backprop backprop_cuda.linkinfo
	cd ./rodinia_3.1/openmp/backprop/opt && rm -rf *.o *~ *.gcda backprop backprop_cuda.linkinfo
	cd ./rodinia_3.1/openmp/lavaMD/org && rm -rf *.o ./kernel/*.o ./kernel/kernel_cpu.gcda ./util/num/*.o ./util/num/*.gcda ./util/timer/*.o ./util/timer/*.gcda ./util/device/*.o ./util/device/*.gcda *.gcda lavaMD
	cd ./rodinia_3.1/openmp/lavaMD/opt && rm -rf *.o ./kernel/*.o ./kernel/kernel_cpu.gcda ./util/num/*.o ./util/num/*.gcda ./util/timer/*.o ./util/timer/*.gcda ./util/device/*.o ./util/device/*.gcda *.gcda lavaMD
	cd ./rodinia_3.1/openmp/hotspot3D/org && rm -rf 3D.gcda 3D *.txt
	cd ./rodinia_3.1/openmp/hotspot3D/opt && rm -rf 3D.gcda 3D *.txt
	cd ./rodinia_3.1/openmp/srad_v2/org && rm -rf srad srad.gcda
	cd ./rodinia_3.1/openmp/srad_v2/opt && rm -rf srad srad.gcda
	cd ./bzip2/org && rm -rf *.o *.gcda libbz2.a bzip2 bzip2recover
	cd ./bzip2/opt && rm -rf *.o *.gcda libbz2.a bzip2 bzip2recover
	cd ./NERSC/smb/org/msgrate && rm -rf *.gcda msgrate.o msgrate	
	cd ./NERSC/smb/opt-msb1/msgrate && rm -rf *.gcda msgrate.o msgrate	
	cd ./NERSC/smb/opt-msb2/msgrate && rm -rf *.gcda msgrate.o msgrate	
	cd ./hoard/original/benchmarks/larson && rm -rf *.gcda larson larson-hoard
	cd ./hoard/optimized/benchmarks/larson && rm -rf *.gcda larson larson-hoard
