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
% each, Phi in this example is gaussian with both positive and negative 
% elements
if(any(runList==001))
    load(filename,'*001');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5;
    opt.thresh=1e-6;
    m=[200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u = 10.^[-1 -2 -3 -4 -5 -6 -7];
    snr=[inf 1e6 1e5 2 5 10 100 1e3 inf];
    for i=1:7
        opt.m=m(i); opt.snr=inf;
        opt=loadLinear(conf,opt);
        for j=3:6
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC001{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC001','-append');

            %opt.continuation=true;
            %opt.alphaStep='FISTA_L1';
            %FISTAC001{i,j}=lasso(conf.Phi,conf.Phit,...
            %    conf.Psi,conf.Psit,conf.y,initSig,opt);
            %save(filename,'FISTAC001','-append');

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

% vary the number of the measurements and get the corresponding good u for 
% each, Use the matrix with positive elements as used in poisson example
if(any(runList==001))
    load(filename,'*001');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5;
    opt.thresh=1e-6;
    m=[200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u = 10.^[-1 -2 -3 -4 -5 -6 -7];
    snr=[inf 1e6 1e5 2 5 10 100 1e3 inf];
    for i=1:7
        opt.m=m(i); opt.snr=inf;
        opt=loadLinear(conf,opt);
        for j=3:6
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC001{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC001','-append');

            %opt.continuation=true;
            %opt.alphaStep='FISTA_L1';
            %FISTAC001{i,j}=lasso(conf.Phi,conf.Phit,...
            %    conf.Psi,conf.Psit,conf.y,initSig,opt);
            %save(filename,'FISTAC001','-append');

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
    load(filename,'*002');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    for j=1:10
        for i=1:7
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.m=m(i); opt.snr=inf; opt.u = u(i);
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            if(j==1) continue; end

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC002{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC002','-append');
            
            opt.continuation=true;
            opt.alphaStep='FISTA_L1';
            FISTAC002{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'FISTAC002','-append');

            opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
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
end

% vary the SNR of measurements, with continuation (continuation is good) to 
% find their corresponding good u, m=600;
if(any(runList==004))
    load(filename,'*004');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5; opt.thresh=1e-6;
    m=[600];
    snr=[1 2 5 10 20 50 100 200 500 1e3 1e4 1e5 1e6];
    u=[10,1,0.1,1e-2,1e-3,1e-4,1e-5,1e-6];

    for i=4:length(snr)
        opt.m=m; opt.snr=snr(i);
        opt=loadLinear(conf,opt);
        for j=1:1
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC004{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC004','-append');
        end
    end
end

% vary the SNR of measurements, with continuation (continuation is good) to 
% find their corresponding good u, m=700;
if(any(runList==014))
    load(filename,'*014');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5; opt.thresh=1e-6;
    m=[700];
    snr=[10,  50, 100, 200,  500,  1e3,  1e4,  1e5,  1e6];
    u=  [ 1, 0.1, 0.1, 0.1, 1e-1, 1e-2, 1e-2, 1e-3, 1e-4]; % this is for m=600

    for i=1:length(snr)
        opt.m=m; opt.snr=snr(i);
        opt=loadLinear(conf,opt);
        for j=1:3
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(i)*10^(j-2);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC014{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC014','-append');
        end
    end
end

% vary the SNR for m=600, u is picked to give the best result
if(any(runList==005))
    load(filename,'*005');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.debugLevel=0;
    opt.maxItr=1e5; opt.thresh=1e-6;
    m=[600];
    snr=[10,  50, 100, 200,  500,  1e3,  1e4,  1e5,  1e6];
    u=  [ 1, 0.1, 0.1, 0.1, 1e-2, 1e-2, 1e-2, 1e-3, 1e-4]; % this is for m=600
    for j=1:10
        for i=1:length(snr)
            opt.m=m; opt.snr=snr(i); opt.u = u(i);
            opt=loadLinear(conf,opt);
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC005{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC005','-append');

            opt.continuation=true;
            opt.alphaStep='FISTA_L1';
            FISTAC005{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'FISTAC005','-append');

            opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
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
end

% vary the SNR for m=700, u is picked to give the best result
if(any(runList==015))
    load(filename,'*015');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.debugLevel=0;
    opt.maxItr=1e5; opt.thresh=1e-6;
    m=[700];
    snr=[10,  50, 100, 200,  500,  1e3,  1e4,  1e5,  1e6];
    u=  [ 1, 0.1, 0.1, 0.1, 1e-1, 1e-2, 1e-2, 1e-3, 1e-4]; % this is for m=600
    for j=1:1
        for i=1:length(snr)
            opt.m=m; opt.snr=snr(i); opt.u = u(i);
            opt=loadLinear(conf,opt);
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC015{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC015','-append');

            opt.continuation=true;
            opt.alphaStep='FISTA_L1';
            FISTAC015{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'FISTAC015','-append');

            opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npg015{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npg015','-append');

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
            out.opt=opt; spiral015{i,j}=out;
            save(filename,'spiral015','-append');

            A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            fpcas015{i,j}=out; fpcas015{i,j}.alpha = conf.Psi(s);
            fpcas015{i,j}.opt = opt; alphaHat=fpcas015{i,j}.alpha;
            fpcas015{i,j}.RMSE=norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha);
            save(filename,'fpcas015','-append');
        end
    end
end


% Poisson example
% vary the number of measurements, with continuation
if(any(runList==006))
    load(filename,'*006');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5; opt.thresh=1e-6;
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

        npg006{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npg006','-append');

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

% The X-ray CT example, test and find the best u for each prjFull
if(any(runList==007))     % FPCAS
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

if(any(runList==008))     % FPCAS
    load(filename,'out009');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360];
    u=10.^[-6 -5 -5 -5 -5 -5];
    subt=[1e-5 1e-6];
    for j=2;
        for i=[5, 1, 2]
            subtolerance=subt(j); conf.prjFull = prjFull(i); opt.u = u(i);
            conf.prjNum = conf.prjFull/2;
            opt=conf.setup(opt);
            opt.maxItr=2e3;
            opt.thresh=1e-12;

            fprintf('%s, i=%d, j=%d\n','SPIRAL-TAP',i,j);
            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
            initSig = maskFunc(conf.FBP(y),opt.mask~=0);
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
            save(filename,'out009','-append');
        end
    end
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


% Plot the figures, or save the data for gnuplot in the paper
if(any(runList==999))
    load(filename); j=1;

    opt.a=0;
    opt=loadLinear(ConfigCT(),opt);
    signal=opt.trueAlpha;
    save('skyline.data','signal','-ascii');

    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200

    npgTime=showResult(npg002,2,'time');
    npgCTime=showResult(npgC002,2,'time');
    fistaCTime=showResult(FISTAC002,2,'time');
    spiralTime=showResult(spiral002,2,'time');
    fpcasTime=showResult(fpcas002,2,'cpu');

    npgCost=showResult(npg002,2,'cost');
    npgCCost=showResult(npgC002,2,'cost');
    fistaCCost=showResult(FISTAC002,2,'cost');
    spiralCost=showResult(spiral002,2,'cost');
    fpcasCost=showResult(fpcas002,2,'f');

    npgRMSE=showResult(npg002,2,'RMSE');
    npgCRMSE=showResult(npgC002,2,'RMSE');
    fistaCRMSE=showResult(FISTAC002,2,'RMSE');
    spiralRMSE=showResult(spiral002,2,'reconerror');
    fpcasRMSE=showResult(fpcas002,2,'RMSE');

    npgTime(:,11:end)=[];
    npgCTime(:,11:end)=[];
    fistaCTime(:,11:end)=[];
    spiralTime(:,11:end)=[];
    fpcasTime(:,11:end)=[];
    npgCost(:,11:end)=[];
    npgCCost(:,11:end)=[];
    fistaCCost(:,11:end)=[];
    spiralCost(:,11:end)=[];
    fpcasCost(:,11:end)=[];
    npgRMSE(:,11:end)=[];
    npgCRMSE(:,11:end)=[];
    fistaCRMSE(:,11:end)=[];
    spiralRMSE(:,11:end)=[];
    fpcasRMSE(:,11:end)=[];

    forSave=[];
    forSave=[forSave, sum(npgTime,2)./sum(npgTime>0,2)];
    forSave=[forSave, sum(npgCTime,2)./sum(npgCTime>0,2)];
    forSave=[forSave, sum(fistaCTime,2)./sum(fistaCTime>0,2)];
    forSave=[forSave, sum(spiralTime,2)./sum(spiralTime>0,2)];
    forSave=[forSave, sum(fpcasTime,2)./sum(fpcasTime>0,2)];

    forSave=[forSave, sum(npgCost,2)./sum(npgCost>0,2)];
    forSave=[forSave, sum(npgCCost,2)./sum(npgCCost>0,2)];
    forSave=[forSave, sum(fistaCCost,2)./sum(fistaCCost>0,2)];
    forSave=[forSave, sum(spiralCost,2)./sum(spiralCost>0,2)];
    forSave=[forSave, sum(fpcasCost,2)./sum(fpcasCost>0,2)];

    forSave=[forSave, sum(npgRMSE,2)./sum(npgRMSE>0,2)];
    forSave=[forSave, sum(npgCRMSE,2)./sum(npgCRMSE>0,2)];
    forSave=[forSave, sum(fistaCRMSE,2)./sum(fistaCRMSE>0,2)];
    forSave=[forSave, sum(spiralRMSE,2)./sum(spiralRMSE>0,2)];
    forSave=[forSave, sum(fpcasRMSE,2)./sum(fpcasRMSE>0,2)];
    forSave=[forSave, m(:)];
    save('varyMeasurement.data','forSave','-ascii');

    clear *Time *Cost *RMSE forSave
    snr=[10 50 100 200 500 1e3 1e4 1e5 1e6];

    npgTime=showResult(npg005,2,'time');
    npgCTime=showResult(npgC005,2,'time');
    fistaCTime=showResult(FISTAC005,2,'time');
    spiralTime=showResult(spiral005,2,'time');
    fpcasTime=showResult(fpcas005,2,'cpu');

    npgCost=showResult(npg005,2,'cost');
    npgCCost=showResult(npgC005,2,'cost');
    fistaCCost=showResult(FISTAC005,2,'cost');
    spiralCost=showResult(spiral005,2,'cost');
    fpcasCost=showResult(fpcas005,2,'f');

    npgRMSE=showResult(npg005,2,'RMSE');
    npgCRMSE=showResult(npgC005,2,'RMSE');
    fistaCRMSE=showResult(FISTAC005,2,'RMSE');
    spiralRMSE=showResult(spiral005,2,'reconerror');
    fpcasRMSE=showResult(fpcas005,2,'RMSE');

    forSave=[];
    forSave=[forSave, sum(npgTime,2)./sum(npgTime>0,2)];
    forSave=[forSave, sum(npgCTime,2)./sum(npgCTime>0,2)];
    forSave=[forSave, sum(fistaCTime,2)./sum(fistaCTime>0,2)];
    forSave=[forSave, sum(spiralTime,2)./sum(spiralTime>0,2)];
    forSave=[forSave, sum(fpcasTime,2)./sum(fpcasTime>0,2)];

    forSave=[forSave, sum(npgCost,2)./sum(npgCost>0,2)];
    forSave=[forSave, sum(npgCCost,2)./sum(npgCCost>0,2)];
    forSave=[forSave, sum(fistaCCost,2)./sum(fistaCCost>0,2)];
    forSave=[forSave, sum(spiralCost,2)./sum(spiralCost>0,2)];
    forSave=[forSave, sum(fpcasCost,2)./sum(fpcasCost>0,2)];

    forSave=[forSave, sum(npgRMSE,2)./sum(npgRMSE>0,2)];
    forSave=[forSave, sum(npgCRMSE,2)./sum(npgCRMSE>0,2)];
    forSave=[forSave, sum(fistaCRMSE,2)./sum(fistaCRMSE>0,2)];
    forSave=[forSave, sum(spiralRMSE,2)./sum(spiralRMSE>0,2)];
    forSave=[forSave, sum(fpcasRMSE,2)./sum(fpcasRMSE>0,2)];
    forSave=[forSave, 10*log10(snr(:))];
    save('varySNR.data','forSave','-ascii');

    fprintf('Poisson example\n');
    forSave=[]; t=0;
    out=ist006{7};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=npg006{7};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=spiral006{7};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.reconerror),t)=out.reconerror;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    mincost=reshape(forSave(:,[1,4,7]),[],1); 
    mincost=min(mincost(mincost~=0));
    idx=(forSave(:,1)~=0); forSave(idx,1)=forSave(idx,1)-mincost;
    idx=(forSave(:,4)~=0); forSave(idx,4)=forSave(idx,4)-mincost;
    idx=(forSave(:,7)~=0); forSave(idx,7)=forSave(idx,7)-mincost;
    save('cost_itr.data','forSave','-ascii');

    !cp vary*.data cost_itr.data skyline.data ~/research/myPaper/asilomar2014/

    clear *Time *Cost *RMSE forSave
end
end

