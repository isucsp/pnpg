%
% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.3 $ $Date: Thu 03 Sep 2015 09:35:30 AM CDT
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
    mex cpuPrj.c mPrj.c common/kiss_fft.c -DCPU=1
    if(gpuDeviceCount>0)
        mex COMPFLAGS="/TP" gpuPrj.obj mPrj.c  common/kiss_fft.c "-LC:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v7.0\lib\x64" -lcudart -DGPU=1
    end
end
cd(pathstr)

cd(['utils' filesep 'L-BFGS-B-C' filesep 'Matlab'])
if(isunix)
    !make
elseif(ispc)
    mex -largeArrayDims -O -g -UDEBUG -I../src ...
    lbfgsb_wrapper.c ../src/lbfgsb.c ../src/linesearch.c ../src/subalgorithms.c ...
    ../src/print.c ../src/linpack.c ../src/miniCBLAS.c ../src/timer.c
end
cd(pathstr)

clear a pathstr

