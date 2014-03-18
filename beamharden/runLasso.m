function [conf,opt] = runLasso(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.2 $ $Date: Tue 18 Mar 2014 01:23:49 AM CDT
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
    conf.imageName='phantom'; %'castSim'; %'phantom' %'twoMaterials'; 

    opt.alphaStep='FISTA_L1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    opt=conf.setup(opt);
    opt.u=1e-4;
    opt.thresh=1e-14;
    opt.maxItr=1e3;
    opt.debugLevel=6;
    %conf.y=conf.y+randn(size(conf.y))*sqrt(1e-8*(norm(conf.y(:)).^2)/length(conf.y(:)));
    prefix='Lasso';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig = conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    %initSig = out0.alpha;
    out0=lasso(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out0','-append');
end

% This section is used to compare different methods for lasso *with* 
% non-negativity constraints
if(any(runList==1)) % FISTA_NNL1
    [conf, opt] = defaultInit();
    i=1; j=1;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.prjFull = 360/6;
    conf.prjNum = conf.prjFull/2;
    conf.PhiMode='cpuPrj'; %'basic'; %'filtered'; %'weighted'; %
    conf.imageName='phantom'; %'castSim'; %'phantom' %'twoMaterials'; 

    opt=conf.setup(opt);

    opt.u=1e-4;
    opt.debugLevel=1;
    opt.alphaStep='FISTA_ADMM_NNL1';%'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    %conf.y=conf.y+randn(size(conf.y))*sqrt(1e-8*(norm(conf.y(:)).^2)/length(conf.y(:)));
    %opt=conf.loadLasso(opt);
    prefix='Lasso';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig = conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    opt.maxItr=1000;
    %initSig = out1.alpha;
    opt.u=1e-2;
    for i=1:9
        opt.u=opt.u*0.1;
        out2{i}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        initSig=out2{i}.alpha;
        save(filename,'out2','-append');
    end
end

% This section is used to compare different methods for lasso *without* 
% non-negativity constraints
if(any(runList==2))
    [conf, opt] = defaultInit();
    i=1; j=1;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.prjFull = 360/6;
    conf.prjNum = conf.prjFull/2;
    conf.PhiMode='cpuPrj'; %'basic'; %'filtered'; %'weighted'; %
    conf.imageName='phantom'; %'castSim'; %'phantom' %'twoMaterials'; 

    opt=conf.setup(opt);

    opt.u=1e-4;
    opt.debugLevel=6;
    opt.alphaStep='FISTA_L1';%'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    
    %conf.y=conf.y+randn(size(conf.y))*sqrt(1e-8*(norm(conf.y(:)).^2)/length(conf.y(:)));
    %opt=conf.loadLasso(opt);
    prefix='Lasso';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig = conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    %initSig = out0.alpha;

    out2=lasso(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out2','-append');
end

if(any(runList==3)) % SPIRAL-G
    [conf, opt] = defaultInit();
    i=1; j=1;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.prjFull = 360/6;
    conf.prjNum = conf.prjFull/2;
    conf.PhiMode='cpuPrj'; %'basic'; %'filtered'; %'weighted'; %
    conf.imageName='phantom_1'; %'castSim'; %'phantom' %'twoMaterials'; 

    opt=conf.setup(opt);

    opt.u=1e-4;
    opt.debugLevel=1;
    initSig = conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    %initSig = out1.alpha;
    opt.alphaStep='FISTA_ADMM_NNL1';%'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    subtolerance=1e-6;
    out=[];
    [out.alpha, out.p, out.cost, out.reconerror, out.time, ...
        out.solutionpath] = ...
        SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
        'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
        'initialization',initSig,'maxiter',opt.maxItr,...
        'miniter',0,'stopcriterion',3,...
        'tolerance',opt.thresh,'truth',opt.trueAlpha,...
        'subtolerance',subtolerance,'monotone',1,...
        'saveobjective',1,'savereconerror',1,'savecputime',1,...
        'reconerrortype',2,...
        'savesolutionpath',1,'verbose',100);
    out_phantom_1_spiral_monotone=out;
    save(filename,'out_phantom_1_spiral_monotone','-append');
    out_phantom_1_FISTA_ADMM_NNL1_a8=lasso(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out_phantom_1_FISTA_ADMM_NNL1_a8','-append');
end

if(any(runList==1002))
end

end

% Don't change the default values arbitrarily, it will cause the old code 
% unrunnable
function [conf, opt] = defaultInit()
    conf = ConfigCT();
    conf.maskType='CircleMask'; %'cvxHull'; %'TightMask'; %
    conf.imageName='castSim'; %'phantom' %'twoMaterials'; %'realct'; %'pellet'; %
    conf.PhiMode='cpuPrj'; %'filtered'; %'weighted'; %
    conf.PhiModeGen='parPrj'; %'cpuPrj'; %'basic';

    opt.continuation = false;
    % from 0 to 7, the higher, the more information is printed. Set to 0 to turn off.
    opt.debugLevel = 5;
    opt.thresh=1e-10;
    opt.maxItr=2e3;
    opt.nu=0;           % nonzero nu means to use nonnegativity constraints
    opt.u=1e-4;         % nonzero u means using sparse
    %opt.a=-6.5;  % aArray=-6.8:0.2:-6.2;
    opt.showImg=true;
end

