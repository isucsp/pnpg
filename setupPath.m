%
% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.3 $ $Date: Thu 27 Aug 2015 01:49:19 AM CDT
%
% 0.4: add variable cleaning statements
% 0.3: add the current path
% 0.2: implement by dbstack, so that it can be called from any location

[a,~]=dbstack('-completenames');
a=a(1);
[pathstr,~,~]=fileparts(a.file);
addpath([pathstr filesep 'beamharden' filesep 'subfunctions']);
addpath([pathstr filesep 'beamharden']);
addpath([pathstr filesep 'tomography']);
addpath([pathstr filesep 'rwt']);
addpath([pathstr filesep 'npg']);
addpath([pathstr filesep 'bhc']);
addpath([pathstr filesep 'prj']);
addpath([pathstr filesep 'utils']);
addpath([pathstr filesep 'utils' filesep 'L-BFGS-B-C' filesep 'Matlab']);
addpath([pathstr filesep 'irt' filesep 'nufft']);
addpath([pathstr filesep 'irt' filesep 'systems']);
addpath([pathstr filesep 'irt']);
addpath([pathstr filesep 'irt' filesep 'emission']);
addpath([pathstr filesep 'irt' filesep 'transmission']);
addpath([pathstr filesep 'irt' filesep 'fbp']);
addpath([pathstr filesep 'irt' filesep 'data']);
addpath([pathstr filesep 'irt' filesep 'utilities']);
addpath([pathstr filesep 'others']);
addpath([pathstr filesep 'others' filesep 'FPC_AS']);
addpath([pathstr filesep 'others' filesep 'FPC_AS' filesep 'src']);
addpath([pathstr filesep 'others' filesep 'FPC_AS' filesep 'prob_gen']);
addpath([pathstr filesep 'others' filesep 'FPC_AS' filesep 'prob_gen' filesep 'classes']);
addpath([pathstr filesep 'others' filesep 'glmnet_matlab']);
addpath([pathstr filesep 'others' filesep 'fpc' filesep 'solvers']);
addpath([pathstr filesep 'others' filesep 'fpc' filesep 'solvers' filesep 'utilities']);

cd 'prj'
if(isunix)
    !make mCPUPrj mParPrj solveTriDiag
    if(gpuDeviceCount>0)
        !make mGPUPrj 
    end
elseif(ispc)
    mex solveTridiag.c
    mex mParPrj.c parPrj.c
    mex -v CFLAGS="-std=c99" cpuPrj.c mPrj.c common/kiss_fft.c -DCPU=1
    if(gpuDeviceCount>0)
        mex gpuPrj.o common/kiss_fft.c mPrj.c ...
	    -DGPU=1 -L. -L./common -L/usr/local/cuda/lib64 -lcuda -lcudart -lglut -lGL
    end
end
cd(pathstr)

cd ['utils' filesep 'L-BFGS-B-C' filesep 'Matlab']
if(isunix)
    !make
elseif(ispc)
    mex -largeArrayDims -lm -O -g -UDEBUG -I../src ...
    lbfgsb_wrapper.c ../src/lbfgsb.c ../src/linesearch.c ../src/subalgorithms.c ...
    ../src/print.c ../src/linpack.c ../src/miniCBLAS.c ../src/timer.c
end
cd(pathstr)

clear a pathstr

