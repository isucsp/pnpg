function [conf,opt] = runLasso(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.2 $ $Date: Fri 14 Mar 2014 01:48:12 AM CDT
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0 || ~isempty(runList))
    filename = [mfilename '.mat'];
    if( exist(filename,'file') ) load(filename);
    else save(filename,'filename');
    end
end

if(nargin==0) runList = [0];
elseif(isempty(runList))
    [conf, opt] = defaultInit();
    opt = conf.setup(opt);
end

%%%%%%%%%%%%%%%%%%%%%%%%
if(any(runList==0)) % reserved for debug and for the best result
    [conf, opt] = defaultInit();
    i=1; j=1;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.prjFull = 360/6;
    conf.prjNum = conf.prjFull/2;
    conf.PhiMode='cpuPrj'; %'basic'; %'filtered'; %'weighted'; %
    conf.imageName='phantom_1'; %'castSim'; %'phantom' %'twoMaterials'; 

    opt.alphaStep='FISTA_NNL1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    opt.alphaStep='ADMM_L1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    opt.nu=0;
    opt.u=1e-4;
    opt.thresh=1e-9;
    opt.debugLevel=5;
    opt=conf.setup(opt);
    %opt=conf.loadLasso(opt);
    prefix='Lasso';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig = conf.FBP(conf.y);
    initSig = 0*initSig(opt.mask~=0);
    initSig = opt.trueAlpha;
    keyboard
    %initSig = out0.alpha;
    out0=lasso(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out0','-append');
end
end

% Don't change the default values arbitrarily, it will cause the old code 
% unrunnable
function [conf, opt] = defaultInit()
    conf = ConfigCT();
    conf.maskType='CircleMask'; %'cvxHull'; %'TightMask'; %
    conf.imageName='castSim'; %'phantom' %'twoMaterials'; %'realct'; %'pellet'; %
    conf.PhiMode='basic'; %'filtered'; %'weighted'; %
    conf.PhiModeGen='parPrj'; %'cpuPrj'; %'basic';

    opt.continuation = false;
    % from 0 to 7, the higher, the more information is printed. Set to 0 to turn off.
    opt.debugLevel = 5;
    opt.thresh=1e-12;
    opt.maxItr=2e3;
    opt.nu=0;           % nonzero nu means to use nonnegativity constraints
    opt.u=1e-4;         % nonzero u means using sparse
    %opt.a=-6.5;  % aArray=-6.8:0.2:-6.2;
    opt.showImg=true;
end

