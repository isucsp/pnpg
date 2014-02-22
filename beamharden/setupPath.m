%
% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.3 $ $Date: Fri 21 Feb 2014 06:24:33 PM CST
%
% 0.4: add variable cleaning statements
% 0.3: add the current path
% 0.2: implement by dbstack, so that it can be called from any location

[a,~]=dbstack('-completenames');
a=a(1);
[pathstr,~,~]=fileparts(a.file);
addpath([pathstr '/subfunctions']);
addpath([pathstr '/rwt']);
addpath([pathstr '/FPC_AS']);
addpath([pathstr '/FPC_AS/src']);
addpath([pathstr '/FPC_AS/prob_gen']);
addpath([pathstr '/FPC_AS/prob_gen/classes']);
addpath([pathstr '/prj']);
addpath([pathstr '/irt/nufft']);
addpath([pathstr '/irt/systems']);
addpath([pathstr '/irt/utilities']);
addpath([pathstr]);
clear a pathstr;
