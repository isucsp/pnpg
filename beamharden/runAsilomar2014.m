function [conf,opt] = runAsilomar2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0 || ~isempty(runList))
    filename = [mfilename '.mat'];
    if(~exist(filename,'file') ) save(filename,'filename'); end
end

if(nargin==0) runList = [0];
elseif(isempty(runList))
    conf=ConfigCT(); opt = conf.setup(); return;
end

% runList rules, abc
% a:
%       0: linear example
%       1: wrist example

%%%%%%%%%%%%%%%%%%%%%%%%

% vary the number of the measurements and get the corresponding good u for 
% each
if(any(runList==001))
    load(filename,'out001');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=2e3;
    opt.thresh=1e-6;
    m=[200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u = 10.^[-1 -2 -3 -4 -5 -6 -7];
    snr=[inf 1e6 1e5 2 5 10 100 1e3 inf];
    opt.alphaStep='FISTA_ADMM_NNL1';
    opt.continuation=true;
    for i=1:7
        opt.m=m(i); opt.snr=inf;
        opt=loadLinear(conf,opt);
        for j=3:6
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(j);
            initSig = conf.Phit(conf.y)*0;

            out001{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out001','-append');

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            fpcas001{i,j}=out; fpcas001{i,j}.alpha = conf.Psi(s);
            fpcas001{i,j}.opt = opt; alphaHat=fpcas001{i,j}.alpha;
            fpcas001{i,j}.RMSE=norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha);
            save(filename,'fpcas001','-append');
        end
        %initSig=out001{i,j}.alpha;
    end
end

% vary the number of measurements, with continuation
% The ture signal is 376 zeropadded.
if(any(runList==002))
    load(filename,'out002');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=2e3; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    opt.alphaStep='FISTA_ADMM_NNL1';
    opt.padZero=376;
    j=1;
    for i=1:7
        fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
        opt.m=m(i); opt.snr=inf; opt.u = u(i);
        opt=loadLinear(conf,opt);
        initSig = conf.Phit(conf.y)*0;

        opt.continuation=true;
        npgC002{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npgC002','-append');

        opt.continuation=false;
        npg002{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npg002','-append');

        subtolerance=1e-5;
        [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
            SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',0,...
            'savesolutionpath',0,'verbose',100);
        out.opt=opt; spiral002{i,j}=out;
        save(filename,'spiral002','-append');

        A = @(xx) conf.Phi(conf.Psi(xx));
        At = @(yy) conf.Psit(conf.Phit(yy));
        AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
        [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
        fpcas002{i,j}=out; fpcas002{i,j}.alpha = conf.Psi(s);
        fpcas002{i,j}.opt = opt; alphaHat=fpcas002{i,j}.alpha;
        fpcas002{i,j}.RMSE=norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha);
        save(filename,'fpcas002','-append');
    end
end

% vary the number of measurements, with continuation
if(any(runList==003))
    load(filename,'out003');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=2e3; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    opt.alphaStep='FISTA_ADMM_NNL1';
    j=1;
    for i=1:7
        fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
        opt.m=m(i); opt.snr=inf; opt.u = u(i);
        opt=loadLinear(conf,opt);
        initSig = conf.Phit(conf.y)*0;

        opt.continuation=true;
        npgC003{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npgC003','-append');

        opt.continuation=false;
        npg003{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npg003','-append');

        subtolerance=1e-5;
        [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
            SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',0,...
            'savesolutionpath',0,'verbose',100);
        out.opt=opt; spiral003{i,j}=out;
        save(filename,'spiral003','-append');

        A = @(xx) conf.Phi(conf.Psi(xx));
        At = @(yy) conf.Psit(conf.Phit(yy));
        AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
        [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
        fpcas003{i,j}=out; fpcas003{i,j}.alpha = conf.Psi(s);
        fpcas003{i,j}.opt = opt; alphaHat=fpcas003{i,j}.alpha;
        fpcas003{i,j}.RMSE=norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha);
        save(filename,'fpcas003','-append');
    end
end

% vary the SNR of measurements, with continuation (continuation is good) to 
% find their corresponding good u
if(any(runList==004))
    load(filename,'out004');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=2e3; opt.thresh=1e-6;
    opt.continuation=true;
    m=[600];
    snr=[1 2 5 10 20 50 100 200 500 1e3 1e4 1e5 1e6];
    u=[10,1,0.1,1e-2,1e-3,1e-4,1e-5];
    opt.alphaStep='FISTA_ADMM_NNL1';
    j=1;
    for i=1:length(snr)
        opt.m=m; opt.snr=snr(i);
        opt=loadLinear(conf,opt);
        for j=1:3
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(j);
            initSig = conf.Phit(conf.y)*0;
            out004{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out004','-append');
        end
    end
end

% vary the SNR
if(any(runList==005))
    load(filename,'out005');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=2e3; opt.thresh=1e-6;
    m=[600];
    snr=[10 50 100 200 500 1e3 1e4 1e5 1e6];
    u=  [1,0.1,0.1,0.1,0.1,1e-2,1e-2,1e-3,1e-3]; % this is for m=600
    opt.alphaStep='FISTA_ADMM_NNL1';
    j=1;
    for i=1:length(snr)
        opt.m=m; opt.snr=snr(i); opt.u = u(i);
        opt=loadLinear(conf,opt);
        fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
        initSig = conf.Phit(conf.y)*0;

        opt.continuation=true;
        npgC005{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npgC005','-append');

        opt.continuation=true;
        npg005{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npg005','-append');

        subtolerance=1e-5;
        [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
            SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',0,...
            'savesolutionpath',0,'verbose',100);
        out.opt=opt; spiral005{i,j}=out;
        save(filename,'spiral005','-append');

        A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
        AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
        [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
        fpcas005{i,j}=out; fpcas005{i,j}.alpha = conf.Psi(s);
        fpcas005{i,j}.opt = opt; alphaHat=fpcas005{i,j}.alpha;
        fpcas005{i,j}.RMSE=norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha);
        save(filename,'fpcas005','-append');
    end
end

% vary the number of measurements, with continuation
if(any(runList==006))
    load(filename,'out006');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=2e3; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    opt.alphaStep='FISTA_ADMM_NNL1';
    opt.noiseType='poisson';
    j=1;
    for i=7:7
        fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
        opt.m=m(i); opt.snr=inf; opt.u = 1e-2;
        opt=loadLinear(conf,opt);
        initSig = conf.Phit(conf.y)*0+1;
        fprintf('min=%d, max=%d\n',min(conf.y), max(conf.y));

        %npg006{i,j}=lasso(conf.Phi,conf.Phit,...
        %    conf.Psi,conf.Psit,conf.y,initSig,opt);
        %save(filename,'npg006','-append');

        opt.alphaStep='IST_ADMM_NNL1';
        ist006{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'ist006','-append');

        subtolerance=1e-5;
        [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
            SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','poisson',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',0,...
            'savesolutionpath',0,'verbose',100);
        out.opt=opt; spiral006{i,j}=out;
        save(filename,'spiral006','-append');
    end
end


if(any(runList==904))
    load(filename,'out004');
    conf=ConfigCT();
    opt=loadLinear(conf);
    u = 10.^[-1 -2 -3 -4 -5 -6 -7];
    opt.maxItr=2e3;
    opt.thresh=1e-12;
    m=[500, 600];
    for j= 2
        opt=loadLinear(conf,opt,m(j));
        for i=[2:3]
            fprintf('%s, i=%d, j=%d\n','SPIRAL',i,j);
            opt.u = u(i);
            initSig = conf.Phit(conf.y)*0;
            out=[];
            subtolerance=1e-6;
            [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
                SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                'initialization',initSig,'maxiter',opt.maxItr,...
                'miniter',0,'stopcriterion',3,...
                'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                'subtolerance',subtolerance,'monotone',1,...
                'saveobjective',1,'savereconerror',1,'savecputime',1,...
                'reconerrortype',0,...
                'savesolutionpath',0,'verbose',100);
            out.opt=opt; out004{i,j}=out;
            save(filename,'out004','-append');
        end
    end
end

if(any(runList==102))     % FPCAS
    load(filename,'out002');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    aArray=[-10:-2];
    for i=[1, 4, 5, 6]
        for j=3:length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            opt.a = aArray(j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
            opt=conf.setup(opt);

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            mu=10^opt.a*max(abs(At(conf.y)));
            option.x0=conf.FBP(conf.y);
            option.x0 = conf.Psit(option.x0(opt.mask~=0));
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            out002{i,j}=out; out002{i,j}.alpha = conf.Psi(s);
            out002{i,j}.opt = opt;
            out002{i,j}.RMSE=1-(out002{i,j}.alpha'*opt.trueAlpha/norm(out002{i,j}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out002','-append');
        end
    end
end

if(any(runList==903)) %solve by Back Projection
    load(filename,'out003');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    for i=1:length(prjFull)
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup();
        out003{i}.img=conf.FBP(conf.y);
        out003{i}.alpha=out003{i}.img(opt.mask~=0);
        out003{i}.RSE=norm(conf.y-conf.Phi(out003{i}.alpha))/norm(conf.y);
        out003{i}.RMSE=1-(out003{i}.alpha'*opt.trueAlpha/norm(out003{i}.alpha)/norm(opt.trueAlpha))^2;
        save(filename,'out003','-append');
    end
end

if(any(runList==904))     % FPCAS after linearization
    load(filename,'out004');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    aArray=[-10:-2];
    for i=[1, 4, 5, 6]
        for j=3:length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            opt.a = aArray(j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
            opt=conf.setup(opt);

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
            mu=10^opt.a*max(abs(At(y)));
            option.x0=conf.FBP(y);
            option.x0 = conf.Psit(option.x0(opt.mask~=0));
            [s, out] = FPC_AS(length(At(y)),AO,y,mu,[],option);
            out004{i,j}=out; out004{i,j}.alpha = conf.Psi(s);
            out004{i,j}.opt = opt;
            out004{i,j}.RMSE=1-(out004{i,j}.alpha'*opt.trueAlpha/norm(out004{i,j}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out004','-append');
        end
    end
end

%solve by Back Projection after linearization
if(any(runList==905))
    load(filename,'out005');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    for i=1:length(prjFull)
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup();
        y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
        out005{i}.img=conf.FBP(y);
        out005{i}.alpha=out005{i}.img(opt.mask~=0);
        out005{i}.RSE=norm(y-conf.Phi(out005{i}.alpha))/norm(y);
        out005{i}.RMSE=1-(out005{i}.alpha'*opt.trueAlpha/norm(out005{i}.alpha)/norm(opt.trueAlpha))^2;
        save(filename,'out005','-append');
    end
end

if(any(runList==996)) % b0, known Ie,
    load(filename,'out006');
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    opt.spectBasis = 'b0';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out6{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out6','-append');
    end
end

if(any(runList==7)) % b1, known Ie,
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    opt.spectBasis = 'b1';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out7{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out7','-append');
    end
end

if(any(runList==71)) % b1, known Ie,
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    conf.PhiMode = 'cpuPrj';
    opt.spectBasis = 'b1';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out71{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out71','-append');
    end
end

if(any(runList==8)) % reserved for debug and for the best result
    [conf, opt] = defaultInit();
    i=1; j=1;
    opt.spectBasis = 'b1';
    %opt.maxIeSteps = 100;
    %conf.theta = (0:6:179)';
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    out8{1}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out8','-append');

    conf.PhiMode='basic'; %'filtered'; %'weighted'; %'cpuFanPar'; %
    conf.PhiModeGen='basic'; %'cpuFanPar'; %'filtered'; %'weighted'; %
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    out8{2}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out8','-append');
end

% SPIRAL-TAP after linearization
if(any(runList==009))
    load(filename,'out009');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u=10.^[-6 -5 -5 -5 -5 -5];
    for i=2:2
        fprintf('%s, i=%d, j=%d\n','SPIRAL-TAP',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt.u = u(i);
        opt=conf.setup(opt);
        opt.maxItr=2e3;
        opt.thresh=1e-12;

        y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
        initSig = maskFunc(conf.FBP(y),opt.mask~=0);
        subtolerance=1e-5;
        out=[];
        [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
            SPIRALTAP_mod(y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',2,...
            'savesolutionpath',0,'verbose',100);
        out.opt=opt; out009{i,j}=out;
        keyboard
        save(filename,'out009','-append');
    end
end

if(any(runList==010))
    load(filename,'out010');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u=10.^[-1 -2 -3 -4 -5 -6];
    for i=1
        for j=6
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2; opt.u = u(j);
            opt=conf.setup(opt);
            opt.alphaStep='FISTA_ADMM_NNL1';
            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
            initSig = maskFunc(conf.FBP(y),opt.mask~=0);
            out010{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,y,initSig,opt);
            save(filename,'out010','-append');
        end
    end
end

if(any(runList==11)) % dis, single AS step,
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    for i=6:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out11{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out11','-append');
    end
end

% dis, max 20 (default) AS steps, CastSim, FISTA_ADMM_NNL1(default)
% Also compare the continuation and non-continuation
% This section has very strong connection with 001
if(any(runList==012))
    load(filename,'out012');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u  =  10.^[-5  -4   -4   -4   -4   -4];
    opt.continuation=false; opt.maxIeSteps=1;
    for i=1:-6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end

    j=2; opt.continuation=true;
    for i=1:-6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end

    j=3; opt.continuation=true; opt.skipIe=true;
    for i=2:-6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end
    j=4; opt.continuation=false; opt.skipIe=true;
    for i=1:6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end
end

if(any(runList==13)) % b0, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b0';
    intval = 6:-1:1;
    j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out13{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out13','-append');
    end
end

if(any(runList==14)) % b0, max AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b0';
    opt.maxIeSteps = 100;
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out14{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out14','-append');
    end
end

if(any(runList==15)) % b1, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out15{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out15','-append');
    end
end

if(any(runList==151)) % b1, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    conf.PhiMode = 'cpuPrj';
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out151{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out151','-append');
    end
end

if(any(runList==152)) % b1, single AS step,
    [conf, opt] = defaultInit();
    i=1; j=1;
    opt.spectBasis = 'b1';
    opt=conf.setup(opt);
    prefix='BeamHard';
    opt.maxItr = 8000;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig = out15{end}.alpha;
    opt.Ie = out15{end}.Ie;
    out152{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out152','-append');
end

if(any(runList==16)) % b1, max AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    opt.maxIeSteps = 100;
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out16{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out16','-append');
    end
end

% dis, compare the FISTA_ADMM_NNL1 and NCG_PR for both continuation and 
% non-continuation
if(any(runList==021))
    load(filename,'out021');
    conf=ConfigCT();
    conf.prjFull = 360; conf.prjNum = conf.prjFull/2; opt.u =  1e-4;
    opt=conf.setup(opt); initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
    opt.maxIeSteps=1;
    i=1; j=1; opt.alphaStep='FISTA_ADMM_NNL1';
    out021{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out021','-append');

    i=1; j=2; opt.alphaStep='NCG_PR';
    out021{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out021','-append');

    i=1; j=3; opt.alphaStep='IST_ADMM_NNL1';
    out021{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out021','-append');

    opt.skipIe=true;
    i=2; j=1; opt.alphaStep='FISTA_ADMM_NNL1';
    out021{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out021','-append');

    i=2; j=2; opt.alphaStep='NCG_PR';
    out021{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out021','-append');

    i=2; j=3; opt.alphaStep='IST_ADMM_NNL1';
    out021{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out021','-append');
end

if(any(runList==22)) % dis, max AS step,
    [conf, opt] = defaultInit();
    opt.maxIeSteps = 100;
    intval = 6:-1:1;
    i=1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        if(j==1) aaa=5; else aaa=1; end
        for i=aaa:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out22{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out22','-append');
        end
    end
end

if(any(runList==23)) % b0, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b0';
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=2:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out23{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out23','-append');
        end
    end
end

if(any(runList==24)) % b0, max AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b0';
    opt.maxIeSteps = 100;
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        if(j==1) aaa=4
        else aaa=1; end
        for i=aaa:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out24{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out24','-append');
        end
    end
end

if(any(runList==25)) % b1, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        if(j==1) aaa=2; else aaa=1; end
        for i=aaa:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out25{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out25','-append');
        end
    end
end

if(any(runList==26)) % b1, max AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    opt.maxIeSteps = 100;
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        if(j==1) aaa=3; else aaa=1; end
        for i=aaa:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out26{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out26','-append');
        end
    end
end

if(any(runList==30)) % Huber function test for best muHuber
    [conf, opt] = defaultInit();
    j=1;
    opt.spectBasis = 'dis';
    opt=rmfield(opt,'muLustig');
    muHuber=[logspace(-1,-15,8), 0.2];
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    for i=2:length(muHuber)
        opt.muHuber=muHuber(i);
        opt=conf.setup(opt);
        initSig=out30{i-1,j}.alpha;
        opt.Ie=out30{i-1,j}.Ie;
        prefix='Huber function for different mu''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out30{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out30','-append');
    end
end

if(any(runList==31)) % test for different u for lustig
    [conf, opt] = defaultInit();
    j=1;
    opt.spectBasis = 'dis';
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    aArray=[-6.5, -9:-4];
    for i=1:length(aArray)
        opt.a = aArray(i);
        opt=conf.setup(opt);
        prefix='continuous method for different u''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out31{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        initSig=out31{i,j}.alpha;
        opt.Ie=out31{i,j}.Ie;
        save(filename,'out31','-append');
    end
end

if(any(runList==32)) % Huber function test for best muHuber
    [conf, opt] = defaultInit();
    j=1;
    opt.spectBasis = 'dis';
    opt=rmfield(opt,'muLustig');
    muHuber=[1e-16 logspace(-1,-15,8), 1e-16];
    muHuber=muHuber(end:-1:1);
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    for i=1:length(muHuber)
        opt.muHuber=muHuber(i);
        opt=conf.setup(opt);
        prefix='Huber function for different mu''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out32{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        initSig=out32{i,j}.alpha;
        opt.Ie=out32{i,j}.Ie;
        save(filename,'out32','-append');
    end
end

% Plot the figures, or save the data for gnuplot in the paper
if(any(runList==999))
    load(filename); j=1;

    opt.a=0;
    opt=loadLinear(ConfigCT(),opt);
    signal=opt.trueAlpha;
    save('skyline.data','signal','-ascii');

    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    for j=1:size(npg002,2)
        for i=1:size(npg002,1)
            npgTime(i,j)=npg002{i,j}.time(end);
            npgCTime(i,j)=npgC002{i,j}.time(end);
            spiralTime(i,j)=spiral002{i,j}.time(end);
            fpcasTime(i,j)=fpcas002{i,j}.cpu(end);

            npgCost(i,j)=npg002{i,j}.cost(end);
            npgCCost(i,j)=npgC002{i,j}.cost(end);
            spiralCost(i,j)=spiral002{i,j}.cost(end);
            fpcasCost(i,j)=fpcas002{i,j}.f(end);

            npgRMSE(i,j)=npg002{i,j}.RMSE(end);
            npgCRMSE(i,j)=npgC002{i,j}.RMSE(end);
            spiralRMSE(i,j)=spiral002{i,j}.reconerror(end);
            fpcasRMSE(i,j)=fpcas002{i,j}.RMSE(end);
        end
    end
    forSave=[npgTime, npgCTime, spiralTime, fpcasTime,...
        npgCost, npgCCost, spiralCost, fpcasCost, ...
        npgRMSE, npgCRMSE, spiralRMSE, fpcasRMSE, m(:)];
    save('varyMeasurement.data','forSave','-ascii');

    clear *Time *Cost *RMSE forSave
    snr=[10 50 100 200 500 1e3 1e4 1e5 1e6];
    for i=1:size(npg005,1)
        npgTime(i,1)=npg005{i,j}.time(end);
        npgCTime(i,1)=npgC005{i,j}.time(end);
        spiralTime(i,1)=spiral005{i,j}.time(end);
        fpcasTime(i,1)=fpcas005{i,j}.cpu(end);

        npgCost(i,1)=npg005{i,j}.cost(end);
        npgCCost(i,1)=npgC005{i,j}.cost(end);
        spiralCost(i,1)=spiral005{i,j}.cost(end);
        fpcasCost(i,1)=fpcas005{i,j}.f(end);

        npgRMSE(i,1)=npg005{i,j}.RMSE(end);
        npgCRMSE(i,1)=npgC005{i,j}.RMSE(end);
        spiralRMSE(i,1)=spiral005{i,j}.reconerror(end);
        fpcasRMSE(i,1)=fpcas005{i,j}.RMSE(end);
    end
    forSave=[npgTime, npgCTime, spiralTime, fpcasTime,...
        npgCost, npgCCost, spiralCost, fpcasCost, ...
        npgRMSE, npgCRMSE, spiralRMSE, fpcasRMSE, 10*log(snr(:))];
    save('varySNR.data','forSave','-ascii');

    !cp vary*.data skyline.data ~/research/myPaper/asilomar2014/

    clear *Time *Cost *RMSE forSave
    
end

end

