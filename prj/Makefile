
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
#CFLAGS = -Wall -c -O2 -fPIC $(INC)
CFLAGS = -std=c99 -fPIC -D_GNU_SOURCE -pthread -fexceptions -fno-omit-frame-pointer
#LFLAGS = -Wall
#COPTIMFLAGS = -O3 -funroll-loops -msse2
 
all: mParPrj solveTriDiag mCPUPrj mGPUPrj cpuPrj

mGPUPrj: mPrj.c gpuPrj.o prj.h
	$(MEX) gpuPrj.o mPrj.c -DGPU=1 -L. -L/usr/local/cuda/lib64 -lcuda -lcudart -lglut -lGL

mCPUPrj: mPrj.c cpuPrj.c common/kiss_fft.c prj.h
	$(MEX) CFLAGS="$(CFLAGS)" cpuPrj.c mPrj.c common/kiss_fft.c -DCPU=1 # -lglut -lGL

mParPrj: mParPrj.c parPrj.c parPrj.h
	$(MEX) CFLAGS="$(CFLAGS)" mParPrj.c parPrj.c

cpuPrj: cpuPrj.c common/kiss_fft.c prj.h
	g++ -Wall -fPIC cpuPrj.c common/kiss_fft.c -o cpuPrj -pthread # -lglut -lGL

gpuPrj.o: gpuPrj.cu prj.h
	nvcc -Xcompiler -fPIC -c -arch=sm_21 gpuPrj.cu

gpuPrj: gpuPrj.o
	nvcc -lglut -lGL gpuPrj.o -o gpuPrj

solveTriDiag:
	$(MEX) solveTridiag.c

clean:
	-rm -vrf *.o *.mex* *~ *.data *.out profiler cpuPrj gpuPrj

