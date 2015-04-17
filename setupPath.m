%
% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.3 $ $Date: Fri 17 Apr 2015 04:55:04 PM CDT
%
% 0.4: add variable cleaning statements
% 0.3: add the current path
% 0.2: implement by dbstack, so that it can be called from any location

[a,~]=dbstack('-completenames');
a=a(1);
[pathstr,~,~]=fileparts(a.file);
addpath([pathstr '/beamharden/subfunctions']);
addpath([pathstr '/beamharden']);
addpath([pathstr '/tomography']);
addpath([pathstr '/rwt']);
addpath([pathstr '/npg']);
addpath([pathstr '/bhc']);
addpath([pathstr '/prj']);
addpath([pathstr '/utils']);
addpath([pathstr '/irt/nufft']);
addpath([pathstr '/irt/systems']);
addpath([pathstr '/irt']);
addpath([pathstr '/irt/emission']);
addpath([pathstr '/irt/transmission']);
addpath([pathstr '/irt/fbp']);
addpath([pathstr '/irt/data']);
addpath([pathstr '/irt/utilities']);
addpath([pathstr '/others/']);
addpath([pathstr '/others/FPC_AS']);
addpath([pathstr '/others/FPC_AS/src']);
addpath([pathstr '/others/FPC_AS/prob_gen']);
addpath([pathstr '/others/FPC_AS/prob_gen/classes']);
addpath([pathstr '/others/glmnet_matlab/']);
addpath([pathstr '/others/fpc/solvers']);
addpath([pathstr '/others/fpc/solvers/utilities']);

cd 'prj'
if(isunix)
    !make mCPUPrj mParPrj solveTriDiag
    if(gpuDeviceCount>0)
        !make mGPUPrj 
    end
end
cd(pathstr)

clear a pathstr hasgpu

