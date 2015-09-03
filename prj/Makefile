
CUDA  := /usr/local/cuda
NVCC  := $(CUDA)/bin/nvcc
MATLAB := $(shell matlab -e | sed -n 's/MATLAB=//p')
MEX := $(MATLAB)/bin/mex
EXT := $(shell $(MATLAB)/bin/mexext)
 
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
 
all: cpu gpu

cpu: mCPUPrj cpuPrj solveTriDiag parPrj.$(EXT)

gpu: mGPUPrj gpuPrj


# compiling gpuPrj

mGPUPrj: $(addsuffix .$(EXT),gpuPrj)

$(addsuffix .$(EXT),gpuPrj) : %.$(EXT) : mPrj.c gpuPrj.o common/kiss_fft.c
	$(MEX) CFLAGS="$(CFLAGS)" -output $* $^ \
	    -DGPU=1 -L. -L./common -L/usr/local/cuda/lib64 -lcudart -lglut -lGL

gpuPrj.o: gpuPrj.cu prj.h utils.h common/kiss_fft.c common/cpu_bitmap.h common/thread.h
	$(NVCC) -Xcompiler -fPIC -c -arch=sm_21 $(filter-out %.h,$^)

gpuPrj: gpuPrj.o
	$(NVCC) -lglut -lGL $< common/kiss_fft.c -o $@
	
# compiling cpuPrj 

mCPUPrj: $(addsuffix .$(EXT),cpuPrj)

$(addsuffix .$(EXT),cpuPrj) : %.$(EXT) : mPrj.c cpuPrj.c common/kiss_fft.c prj.h utils.h common/thread.h
	$(MEX) CFLAGS="$(CFLAGS)" -output $* $(filter-out %.h,$^) -DCPU=1 # -lglut -lGL

cpuPrj: cpuPrj.c common/kiss_fft.c prj.h utils.h common/cpu_bitmap.h common/thread.h
	g++ -Wall -fPIC $(filter-out %.h,$^) -o $@ -pthread -lglut -lGL

# compiling parallel projection operator 

$(addsuffix .$(EXT),parPrj) : %.$(EXT) : mParPrj.c parPrj.c parPrj.h common/thread.h
	$(MEX) CFLAGS="$(CFLAGS)" -output $* $(filter-out %.h,$^)

solveTriDiag: $(solveTriDiag .$(EXT),solveTriDiag)

$(solveTriDiag .$(EXT),solveTriDiag) : solveTriDiag.c
	$(MEX) $^

.PHONY : clean
clean:
	-rm -vrf *.o *.mex* *~ *.data *.out profiler cpuPrj gpuPrj

