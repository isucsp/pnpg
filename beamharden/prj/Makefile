
CUDA  := /usr/local/cuda
MATLAB = $(shell matlab -e | sed -n 's/MATLAB=//p')

MEX = $(MATLAB)/bin/mex
CXX = g++
AR = ar
 
INC_CU = -I$(CUDA)/include -I$(CUDA)/samples/common/inc
INC_MATLAB = -I$(MATLAB)/extern/include
LIBS = -L$(CUDA)/lib64 -L$(CUDA)/lib -Wl,-rpath,$(CUDA)/lib64
LIBS_CU = -lcudart -lcusparse -lcublas -lcufft -lcurand
LIBS_GL = -lglut -lGL
 
# Regular C++ part
CFLAGS = -Wall -c -O2 -fPIC $(INC)
CFLAGS = -fPIC -D_GNU_SOURCE -pthread -fexceptions -arch glnxa64 
LFLAGS = -Wall
COPTIMFLAGS = -O3 -funroll-loops -msse2
 
all: mParPrj solveTriDiag mCPUPrj cpuPrj

mGPUPrj: mPrj.c gpuPrj.o
	$(MEX) mPrj.c gpuPrj.o -DGPU=1 -L. -L/usr/local/cuda/lib64 -lcuda -lcudart -lglut -lGL

mCPUPrj: mPrj.c cpuPrj.o
	$(MEX) mPrj.c cpuPrj.o -DCPU=1

mParPrj: mParPrj.c parPrj.c parPrj.h
	$(MEX) mParPrj.c parPrj.c

cpuPrj.o: cpuPrj.c prj.h
	g++ -Wall -fPIC -c cpuPrj.c -o cpuPrj.o

cpuPrj: cpuPrj.c common/kiss_fft.c
	g++ cpuPrj.c common/kiss_fft.c -o cpuPrj -pthread -lglut -lGL

gpuPrj.o: gpuPrj.cu prj.h
	nvcc -Xcompiler -fPIC -c -arch=sm_21 gpuPrj.cu

gpuPrj: gpuPrj.o
	nvcc -lglut -lGL gpuPrj.o -o gpuPrj

solveTriDiag:
	$(MEX) solveTridiag.c

clean:
	-rm -vrf *.o *.mex* *~ *.data *.out

