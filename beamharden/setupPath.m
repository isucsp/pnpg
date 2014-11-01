%
% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.3 $ $Date: Fri 24 Oct 2014 01:49:04 AM CDT
%
% 0.4: add variable cleaning statements
% 0.3: add the current path
% 0.2: implement by dbstack, so that it can be called from any location

[a,~]=dbstack('-completenames');
a=a(1);
[pathstr,~,~]=fileparts(a.file);
addpath([pathstr '/subfunctions']);
addpath([pathstr '/../rwt']);
addpath([pathstr '/../npg']);
addpath([pathstr '/../prj']);
addpath([pathstr '/../irt/nufft']);
addpath([pathstr '/../irt/systems']);
addpath([pathstr '/../irt/utilities']);
addpath([pathstr '/../others/']);
addpath([pathstr '/../others/FPC_AS']);
addpath([pathstr '/../others/FPC_AS/src']);
addpath([pathstr '/../others/FPC_AS/prob_gen']);
addpath([pathstr '/../others/FPC_AS/prob_gen/classes']);
addpath([pathstr '/../others/glmnet_matlab/']);
addpath([pathstr '/../others/fpc/solvers']);
addpath([pathstr '/../others/fpc/solvers/utilities']);
addpath([pathstr]);

% cd '../prj'
% %!make clean
% !make mCPUPrj mParPrj solveTriDiag mGPUPrj 
% cd(pathstr)
% clear a pathstr;

