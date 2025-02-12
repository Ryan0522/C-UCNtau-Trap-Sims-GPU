# Variables
OPT=-O3 -pg -g
CPP=mpic++
CPPFLAGS=$(OPT) -std=c++11
CC=mpicc
CFLAGS=$(OPT)

NVCC=nvcc
NVCCFLAGS=$(OPT)

# Check for MPI compilers
MPICPP := $(shell command -v mpic++ 2> /dev/null)
MPICC := $(shell command -v mpicc 2> /dev/null)
ifndef MPICPP
	CPP=CC
endif
ifndef MPICC
	CC=cc
endif

# Files: assume your GPU files now have .cu extension
# CPU/C++ files:
CPU_OBJS= bin/xorshift.o bin/quant_refl.o bin/geometry.o bin/track_gen.o
# CUDA files:
CUDA_OBJS= bin/fields_nate.o bin/symplectic.o bin/trackUCN.o bin/sim.o

all: sim

debug: OPT=-g
debug: all

# Compile CPU source files with mpicc or mpic++
bin/xorshift.o: src/xorshift.c inc/xorshift.h
	$(CC) $(CFLAGS) -c -o $@ $<

bin/quant_refl.o: src/quant_refl.cpp inc/quant_refl.hpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

bin/geometry.o: src/geometry.cpp inc/geometry.hpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

bin/track_gen.o: src/track_gen.cpp inc/track_gen.hpp
	$(CPP) $(CPPFLAGS) -c -o $@ $<

# Compile CUDA source files with nvcc.
bin/fields_nate.o: src/fields_nate.cu inc/fields_nate.h
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

bin/symplectic.o: src/symplectic.cu inc/symplectic.hpp inc/fields_nate.h
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

bin/trackUCN.o: src/trackUCN.cu inc/trackUCN.hpp
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

bin/sim.o: sim.cu
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<

# Link all objects (using mpic++ for MPI linking)
sim: $(CPU_OBJS) $(CUDA_OBJS)
	$(CPP) $(LFLAGS) -pg -g -o ./in/PENdebug1550 $(CPU_OBJS) $(CUDA_OBJS)

# Optional rule for find_min (if needed)
bin/find_min.o: find_min.c
	$(CC) $(CFLAGS) -c -o $@ $<

find_min: bin/find_min.o bin/fields_nate.o
	$(CC) $(CFLAGS) $^ -o find_min $(LFLAGS) -lgsl

clean:
	find ./bin/ -type f -name '*.o' -delete
	rm -rf sim find_min ./in/PENdebug1550
