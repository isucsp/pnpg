
CUDA  := /usr/local/cuda
MATLAB = /usr/local/MATLAB/R2013b

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
 
all: mParPrj solveTriDiag mPrj gpuPrj

mPrj: mPrj.c gpuPrj.o
	$(MEX) mPrj.c gpuPrj.o -L. -L/usr/local/cuda/lib64 -lcuda -lcudart -lglut -lGL

#$(MEX) mPrj.c gpuPrj.o -L. -L/usr/local/cuda/lib64 -lcuda -lcudart -lglut -lGL 

mParPrj: mParPrj.c parPrj.o
	$(MEX) mParPrj.c parPrj.o 

parPrj.o: parPrj.c parPrj.h
	gcc -c -fPIC -lpthread parPrj.c

cpuPrj: cpuPrj.o
	g++ cpuPrj.o -o cpuPrj -lglut -pthread

cpuPrj.o: cpuPrj.c prj.h
	g++ -Wall -c cpuPrj.c -o cpuPrj.o

gpuPrj: gpuPrj.o
	nvcc -lglut -lGL gpuPrj.o -o gpuPrj

gpuPrj.o: gpuPrj.cu prj.h 
	nvcc -Xcompiler -fPIC -c -arch=sm_21 gpuPrj.cu

solveTriDiag:
	$(MEX) solveTridiag.c

gpuRun: gpuPrj.o
	nvcc -lglut -lGL gpuPrj.o -o gpuPrj -run

clean:
	-rm -vrf *.o *.mex* *~ *.data

