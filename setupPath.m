%
% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.3 $ $Date: Tue, Aug 25, 2015  4:12:28 PM
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

end
cd(pathstr)

cd 'utils' filesep 'L-BFGS-B-C' filesep 'Matlab'
if(isunix)
    !make
elseif(ispc)

end
cd(pathstr)

clear a pathstr hasgpu

