function [conf,opt] = runAsilomar2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0) runList = [0];
elseif(isempty(runList))
    conf=ConfigCT(); opt = conf.setup(); return;
end

% runList rules, abc
% a:
%       0: linear example
%       1: wrist example

%%%%%%%%%%%%%%%%%%%%%%%%

if(any(runList==000))
    filename = [mfilename '_000.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    n=1000;p=10000;
    nzc=fix(p/10);
    x=randn(n,p);
    iz=randsample(n*p,n*p*0.85,false);
    x(iz)=0;
    sx=sparse(x);
    beta=randn(nzc,1);
    fx=x(:,1:nzc)*beta;
    eps=randn(n,1);
    y=fx+eps;
    tic;
    gns=glmnet(sx,y,'gaussian',options);
    gns.time=toc;
    conf.Phi=@(aaa) x*aaa;
    conf.Phit=@(aaa) x'*aaa;
    conf.Psi=@(aaa) aaa;
    conf.Psit=@(aaa) aaa;
    opt.trueAlpha=[beta(:);zeros(p-nzc,1)];
    opt.thresh=1e-12;
    for i=1:length(gn.lambda)
        opt.u=gn.lambda(i)*n;
        gn.rmse(i)=sqrNorm(gn.beta(:,i)-opt.trueAlpha)/trueAlphaNorm;
        gn.cost(i)=0.5*sqrNorm(y-gn.a0(i)-x*gn.beta(:,i))+opt.u*pNorm(gn.beta(:,i),1);
        npgsc{i}=NPG.NPGsc(conf.Phi,conf.Phit,conf.Psi,conf.Psit,y-gn.a0(i),zeros(p,1),opt);
        npgs{i}=NPG.NPGs(conf.Phi,conf.Phit,conf.Psi,conf.Psit,y-gn.a0(i),zeros(p,1),opt);
        fpcas{i}=NPG.FPCas(conf.Phi,conf.Phit,conf.Psi,conf.Psit,y-gn.a0(i),zeros(p,1),opt);
    end
    save(filename);
end

if(any(runList==001))
    filename = [mfilename '_001.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT(); k=1;
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    for i=1:length(m)
        opt.m=m(i); opt.snr=inf;
        opt=loadLinear(conf,opt);
        initSig = conf.Phit(conf.y)*0;
        u = (1e-5)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
        for j=1:4
            fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
            opt.u = u*10^(j-2);

            opt.continuation=false; opt.alphaStep='FISTA_ADMM_NNL1';
            npg{i,j,k}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npg','-append');

            opt.continuation=false; opt.alphaStep='FISTA_L1';
            npgs{i,j,k}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgs','-append');

            ppsi = @(yyy,uuu) conf.Psi(Utils.softThresh(conf.Psit(yyy),uuu));
            rrrr = @(xxx) pNorm(conf.Psit(xxx),1);
            [x_SpaRSA,x_debias_SpaRSA,obj_SpaRSA,times_SpaRSA,debias_start_SpaRSA,mse]=...
                SpaRSA_mod(conf.y,conf.Phi,opt.u,...
                'AT',conf.Phit,...
                'Psi',ppsi,...
                'Phi',rrrr,...
                'Safeguard',1,...
                'Initialization',initSig,...
                'StopCriterion',5,...
                'ToleranceA',opt.thresh, ...
                'True_x',opt.trueAlpha,...
                'MaxiterA',opt.maxItr);
            clear('out');
            out.alpha=x_SpaRSA; out.cost=obj_SpaRSA; out.time=times_SpaRSA;
            out.RMSE=mse.mses/sqrNorm(opt.trueAlpha)*length(opt.trueAlpha);
            out.stepSize=mse.stepSize; out.difAlpha=mse.difAlpha;
            sparsa{i,j,k}=out;
            save(filename,'sparsa','-append');
        end
    end
end

% vary the inter loop criteria
if(any(runList==021))
    load(filename,'*021');
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-10; opt.debugLevel=1;
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    i=4;
    opt.m=m(i); opt.snr=inf; opt.u = u(i);
    opt=loadLinear(conf,opt);
    for j=9:2:11
        for i=1
            opt.admmAbsTol=10^(-j);
            opt.admmTol=0;
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npg021{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npg021','-append');
        end
    end
end
if(any(runList==921))
    load(filename,'*021');
    figure(1); t=0;
    style={'r','g','b','k','r-.','g-.','b-.','k-.'};
    cost=showResult(npg021,2,'cost');
    cost=min(cost(cost>0));
    for j=1:2:11
        t=t+1;
        subplot(4,2,1); semilogy(npg021{1,j}.difAlpha,style{t}); hold on;
        subplot(4,2,3); semilogy(npg021{1,j}.RMSE,style{t}); hold on;
        subplot(4,2,5); semilogy(npg021{1,j}.time,style{t}); hold on;
        subplot(4,2,2); semilogy(npg021{1,j}.time,npg021{1,j}.difAlpha,style{t}); hold on;
        subplot(4,2,4); semilogy(npg021{1,j}.time,npg021{1,j}.RMSE,style{t}); hold on;
        subplot(4,2,6); plot    (npg021{1,j}.time,npg021{1,j}.time,style{t}); hold on;
        subplot(4,2,7); semilogy(npg021{1,j}.cost-cost,style{t}); hold on;
        subplot(4,2,8); semilogy(npg021{1,j}.time,npg021{1,j}.cost-cost,style{t}); hold on;
    end
    legend('10^{-1}','10^{-3}','10^{-5}','10^{-7}','10^{-9}','10^{-11}');
    subplot(4,2,1); ylabel('\delta x'); ylim([1e-10,1]);
    subplot(4,2,3); ylabel('RSE');
    subplot(4,2,5); ylabel('time');
    subplot(4,2,7); ylabel('cost'); ylim([1e-10,1e2]);
    subplot(4,2,2); xlim([0,200]); ylim([1e-10,1]);
    subplot(4,2,4); xlim([0,200]);
    subplot(4,2,6); xlim([0,200]);
    subplot(4,2,8); xlim([0,200]);
    subplot(4,2,7); xlabel('# iterations');
    subplot(4,2,8); xlabel('time'); ylim([1e-10,1e2]);
end

% vary the inter loop criteria
if(any(runList==031))
    load(filename,'*031');
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-10; opt.debugLevel=1;
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    i=4;
    opt.m=m(i); opt.snr=inf; opt.u = u(i);
    opt=loadLinear(conf,opt);
    for j=1:5
        for i=1
            opt.admmAbsTol=0;
            opt.admmTol=10^(-j);
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npg031{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npg031','-append');
        end
    end
end

if(any(runList==931))
    load(filename,'*031');
    figure; t=0;
    style={'r','g','b','k','r-.','g-.','b-.','k-.'};
    cost=min(reshape(showResult(npg031,2,'cost'),[],1));
    for j=1:5
        t=t+1;
        subplot(4,2,1); semilogy(npg031{1,j}.difAlpha,style{t}); hold on;
        subplot(4,2,3); semilogy(npg031{1,j}.RMSE,style{t}); hold on;
        subplot(4,2,5); semilogy(npg031{1,j}.time,style{t}); hold on;
        subplot(4,2,2); semilogy(npg031{1,j}.time,npg031{1,j}.difAlpha,style{t}); hold on;
        subplot(4,2,4); semilogy(npg031{1,j}.time,npg031{1,j}.RMSE,style{t}); hold on;
        subplot(4,2,6); plot    (npg031{1,j}.time,npg031{1,j}.time,style{t}); hold on;
        subplot(4,2,7); semilogy(npg031{1,j}.cost-cost,style{t}); hold on;
        subplot(4,2,8); semilogy(npg031{1,j}.time,npg031{1,j}.cost-cost,style{t}); hold on;
    end
    legend('10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}');
    subplot(4,2,1); ylabel('\delta x');
    subplot(4,2,3); ylabel('RSE');
    subplot(4,2,5); ylabel('time');
    subplot(4,2,7); ylabel('cost');
    subplot(4,2,2); xlim([0,400]);
    subplot(4,2,4); xlim([0,400]);
    subplot(4,2,6); xlim([0,400]);
    subplot(4,2,8); xlim([0,400]);
    subplot(4,2,7); xlabel('# iterations');
    subplot(4,2,8); xlabel('time');
end

% vary the number of measurements, with continuation
if(any(runList==002))
    filename = [mfilename '_002.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
    m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u = [1e-3,1e-3,1e-4,1e-4,1e-5,1e-5,1e-6,1e-6,1e-6];
    for k=1:5
        for i=1:length(m)
            opt.m=m(i); opt.snr=inf;
            [opt,~,invEAAt,Phi]=loadLinear(conf,opt);
            initSig = conf.Phit(invEAAt*conf.y);

            if(k==2) save(filename); return; end;
            opt.u = u(i)*10.^(-2:2);
            %gnet{i,k}=NPG.glmnet(Phi,wvltMat(length(opt.trueAlpha),conf.dwt_L,conf.daub),conf.y,initSig,opt);

            for j=1:5
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
                opt.u = u(i)*10^(j-3)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);

                % if(i~=6) continue; end
                % if(j>3) continue; end

                % pgc   {i,j,k}=NPG.PGc(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                % continue;

                % temp=opt.thresh; opt.thresh=1e-8;
                % %pgc12{i,j,k}=NPG.PGc(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                % %sparsn12{i,j,k}=NPG.SpaRSAp(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                % opt.thresh=temp;
                % continue;

                % if(k==1 && i==2)
                %     npgsT {i,j,k}=NPG.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

                %     continue;

                %     opt.initStep='fixed';
                %     fistal{i,j,k}=NPG.FISTA(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                %     opt=rmfield(opt,'initStep');
                % end

                % opt.contShrnk=0.25;
                % opt.contCrtrn=1e-3;
                % opt.adaptiveStep=false;

                npgsc    {i,j,k}=NPG.NPGsc    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgc     {i,j,k}=NPG.NPGc     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                continue;

                fpc      {i,j,k}=NPG.FPC      (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgs     {i,j,k}=NPG.NPGs     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgc_nads{i,j,k}=NPG.NPGc_nads(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npg_nads {i,j,k}=NPG.NPG_nads (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                pgc      {i,j,k}=NPG.PGc      (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                sparsa   {i,j,k}=NPG.SpaRSA   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                sparsn   {i,j,k}=NPG.SpaRSAp  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npg      {i,j,k}=NPG.main     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                spiral   {i,j,k}=NPG.SPIRAL   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                fista    {i,j,k}=NPG.FISTA    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                fpcas    {i,j,k}=NPG.FPCas    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
            end

            save(filename);
        end
    end
end

% Poisson identity link example
% vary the number of measurements, with continuation
if(any(runList==004))
    filename = [mfilename '_004.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    conf.imageName = 'wrist';
    conf.PhiMode = 'gpuPrj';
    conf.PhiModeGen = 'gpuPrj';
    conf.dist = 4000;
    conf.beamharden = false;

    prjFull = [60, 80, 100, 120, 180, 360];
    a=        [-2, -2,  -2,  -2,  -2,  -2];

    opt.maxItr=1e4; opt.thresh=1e-6; opt.snr=1e4; opt.debugLevel=1;
    opt.noiseType='poisson';

    j=1;
    for i=1:length(prjFull)
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull;
        opt=conf.setup(opt);
        L = @(aaa) Utils.poissonModel(aaa,conf.Phi,conf.Phit,conf.y);
        [~,g]=L(initSig*0); u_max=pNorm(conf.Psit(g),inf);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

        fprintf('min=%d, max=%d\n',min(conf.y), max(conf.y));

        opt.fullcont=true;
        opt.u=10.^(-2:-0.5:-6)*u_max;
        if(any([1 3 6]==i))
            npgFull   {i}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
            npgsFull  {i}=NPG.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        end
        opt.fullcont=false;

        save(filename);
        continue;

        tic
        fbp{i,1}.img=conf.FBP(conf.y);
        fbp{i,1}.time=toc;
        fbp{i,1}.alpha=fbp{i,j}.img(opt.mask~=0);
        fbp{i,1}.RSE=sqrNorm(conf.y-conf.Phi(fbp{i,j}.alpha))/sqrNorm(conf.y);
        fbp{i,1}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);

        opt.u = u(i)*10^(j-3);
        fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
        npg   {i,j,k}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npgc  {i,j,k}=NPG.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npgs  {i,j,k}=NPG.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        spiral{i,j,k}=NPG.SPIRAL (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

        save(filename);
    end
end


% vary the SNR for m=600, u is picked to give the best result
K = 5;
if(any(runList==005))
    filename = [mfilename '_005.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.debugLevel=0;
    opt.maxItr=1e4; opt.thresh=1e-6;
    m=[600];
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7];
    a  =[1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3,1e-4,1e-4,1e-5];

    for k=1:K
        for i=1:length(snr)
            opt.m=m; opt.snr=snr(i);
            [opt,~,invEAAt]=loadLinear(conf,opt);
            initSig = conf.Phit(invEAAt*conf.y);
            for j=1:5
                opt.u = a(i)*10^(j-3);
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);

                npg   {i,j,k}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgc  {i,j,k}=NPG.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                fista {i,j,k}=NPG.FISTA  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgs  {i,j,k}=NPG.NPGsc  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                spiral{i,j,k}=NPG.SPIRAL (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                fpcas {i,j,k}=NPG.FPCas  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                sparsa{i,j,k}=NPG.SpaRSA (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                sparsn{i,j,k}=NPG.SpaRSAp(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
            end
            save(filename);
        end
    end
end

% Poisson example
% vary the number of measurements, with continuation
if(any(runList==006))
    filename = [mfilename '_006.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2];
    opt.noiseType='poisson'; opt.matrixType='nonneg';
    for k=1:5
        for i=1:length(m)
            opt.m=m(i); opt.snr=inf;
            [opt,~,invEAAt]=loadLinear(conf,opt);
            initSig=conf.Phit(invEAAt*conf.y);
            fprintf('min=%d, max=%d\n',min(conf.y), max(conf.y));
            if(k<2) continue; end;
            for j=1:5
                opt.u = u(i)*10^(j-3);
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
                npgc  {i,j,k}=NPG.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npg   {i,j,k}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgs  {i,j,k}=NPG.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                spiral{i,j,k}=NPG.SPIRAL (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
            end
            save(filename);
        end
    end
end

% log link Poisson example
% vary the number of measurements, with continuation
if(any(runList==007))
    filename = [mfilename '_007.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800, 900]; % should go from 200
    a=[  -2,  -2,  -2,  -2,  -2,  -2,  -2,  -2];
    opt.noiseType='poissonLogLink'; opt.matrixType='nonneg';
    for k=1:5
        for i=1:length(m)
            opt.m=m(i);
            [opt,~,invEAAt,Phi]=loadLinear(conf,opt);
            temp=conf.y; temp(temp==0)=1;
            initSig=-conf.Phit(invEAAt*log(temp/max(conf.y)));
            fprintf('min=%d, max=%d\n',min(conf.y), max(conf.y));
            
            if(k==2) save(filename); return; end
            if(any([2 3 4 6]==i)) continue; end

            opt.fullcont=true;

                L = @(aaa) Utils.poissonModelLogLink(aaa,conf.Phi,conf.Phit,conf.y);
                [~,g]=L(initSig*0);
                u_max=pNorm(conf.Psit(g),inf);

                opt.u=10.^(0:-0.1:-4)*u_max;
                gnet      {i,k}=NPG.glmnet(Phi,wvltMat(length(opt.trueAlpha),conf.dwt_L,conf.daub),conf.y,initSig,opt);
                npgFull   {i,k}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgsFull  {i,k}=NPG.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                
                opt.noiseType='gaussian';

                    temp=conf.y; temp(temp==0)=1; temp=log(opt.I0./temp); temp=temp.*sqrt(conf.y);
                    wPhi=@(xxx) sqrt(conf.y).*conf.Phi(xxx);
                    wPhit=@(xxx) conf.Phit(sqrt(conf.y).*xxx);

                    L = @(aaa) Utils.linearModel(aaa,wPhi,wPhit,temp);
                    [~,g]=L(initSig*0); u_max=pNorm(conf.Psit(g),inf);
                    opt.u=10.^(-1:-0.1:-3)*u_max;
                    npglwFull {i,k}=NPG.main   (wPhi,wPhit,conf.Psi,conf.Psit,temp,initSig,opt);
                    npgslwFull{i,k}=NPG.NPGs  (wPhi,wPhit,conf.Psi,conf.Psit,temp,initSig,opt);

                    temp=conf.y; temp(temp==0)=1; temp=log(opt.I0./temp);
                    L = @(aaa) Utils.linearModel(aaa,conf.Phi,conf.Phit,temp);
                    [~,g]=L(initSig*0); u_max=pNorm(conf.Psit(g),inf);
                    opt.u=10.^(-1:-0.1:-3)*u_max;
                    npglFull {i,k}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,temp,initSig,opt);
                    npgslFull{i,k}=NPG.NPGs  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,temp,initSig,opt);

                opt.noiseType='poissonLogLink';

            opt.fullcont=false;

            continue;

            for j=1:5
                opt.u = u(i)*10^(j-3)*u_max;
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);

                opt.noiseType='gaussian';
                temp=conf.y; temp(temp==0)=1; temp=log(opt.I0./temp); temp=temp.*sqrt(conf.y);
                wPhi=@(xxx) sqrt(conf.y).*conf.Phi(xxx);
                wPhit=@(xxx) conf.Phit(sqrt(conf.y).*xxx);
                npgclw {i,j,k}=NPG.NPGc   (wPhi,wPhit,conf.Psi,conf.Psit,temp,initSig,opt);
                npgsclw{i,j,k}=NPG.NPGsc  (wPhi,wPhit,conf.Psi,conf.Psit,temp,initSig,opt);

                temp=conf.y; temp(temp==0)=1; temp=log(opt.I0./temp);
                npgcl {i,j,k}=NPG.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,temp,initSig,opt);
                npgscl{i,j,k}=NPG.NPGsc  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,temp,initSig,opt);
                opt.noiseType='poissonLogLink';

                npgc  {i,j,k}=NPG.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npg   {i,j,k}=NPG.main   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgsc {i,j,k}=NPG.NPGsc  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npgs  {i,j,k}=NPG.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

            end
            save(filename);
        end
    end
end

% The X-ray CT example, test and find the best u for each prjFull
if(any(runList==008))     % FPCAS
    filename = [mfilename '_008.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt'); RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    noise = randn(360*1024,1);
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    prjFull = [60, 80, 100, 120, 180, 360];
    opt.debugLevel=1;
    j=1;
    for i=1:length(prjFull)
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        conf.y=conf.Phi(opt.trueAlpha); % equivalent to linear projection
        v = (norm(conf.y)/sqrt(opt.snr*length(conf.y)));
        conf.y=conf.y+v*noise(1:length(conf.y));
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

        % tic
        % fbp{i,j}.img=conf.FBP(conf.y);
        % fbp{i,j}.time=toc;
        % fbp{i,j}.alpha=fbp{i,j}.img(opt.mask~=0);
        % fbp{i,j}.RSE=sqrNorm(conf.y-conf.Phi(fbp{i,j}.alpha))/sqrNorm(conf.y);
        % %fbp{i,j}.RMSE=1-(fbp{i,j}.alpha'*opt.trueAlpha/norm(fbp{i,j}.alpha)/norm(opt.trueAlpha))^2;
        % fbp{i,j}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        % fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);

        fprintf('%s, i=%d, j=%d\n','X-ray CT example',i,j);
        if(j==1)
            u=10.^[-6 -5 -5 -5 -5 -5]; % snr=inf  for the first column
            opt.maxItr=2e3; opt.thresh=1e-6;
            opt.u=u(i);
        end
        if(j==2)
            u=10.^[-8 -8 -8 -8 -8 -9]; % snr=inf  for the first column
            opt.maxItr=5e3; opt.thresh=1e-7;
            opt.u = u(i)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
        end
        if(j==3)
            u=10.^[-9 -9 -9 -9 -9 -10]; % snr=inf for the 3rd column
            opt.maxItr=1e4; opt.thresh=1e-8;
            opt.u = u(i)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
        end

        sparsn   {i,j}=NPG.SpaRSAp  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename);
        continue
        npgc_nads{i,j}=NPG.NPGc_nads(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npg_nads {i,j}=NPG.NPG_nads (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npg      {i,j}=NPG.main     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npgs     {i,j}=NPG.NPGs     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

        npgsc    {i,j}=NPG.NPGsc    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npgc     {i,j}=NPG.NPGc     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

        spiral   {i,j}=NPG.SPIRAL   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        fpcas    {i,j}=NPG.FPCas    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        sparsa   {i,j}=NPG.SpaRSA   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        fista    {i,j}=NPG.FISTA    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

    end
end

% the last example with glassbeads under the poisson model with log link 
% function
if(any(runList==009))
    filename = [mfilename '_009.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    conf=ConfigCT();
    conf.imageName = 'glassBeadsSim';
    conf.PhiMode = 'gpuPrj';
    conf.PhiModeGen = 'gpuPrj';
    conf.dist = 17000;
    conf.beamharden = false;

    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    opt.maxItr=2e3; opt.thresh=1e-6; opt.snr=1e6; opt.debugLevel=1;
    opt.noiseType='poissonLogLink'; %'gaussian'; %
    %for i=1:length(prjFull)
    for i=3:5
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(-log(conf.y)),opt.mask~=0);
        % initSig = opt.trueAlpha;
        for j=1:1
            fprintf('%s, i=%d, j=%d\n','X-ray CT example glassBeads Simulated',i,j);

            % fbp{i,j}.img=conf.FBP(-log(conf.y));
            % fbp{i,j}.alpha=fbp{i,j}.img(opt.mask~=0);
            % fbp{i,j}.RSE=sqrNorm(conf.y-conf.Phi(fbp{i,j}.alpha))/sqrNorm(conf.y);
            % fbp{i,j}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
            % fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);
            % save(filename,'fbp','-append');

            % u=10.^[-7 -7 -7 -6 -5 -5];
            % opt.u = u(i)*10^(j-2);
            % opt.continuation = false; opt.alphaStep='FISTA_ADMM_NNL1';
            % npg{i,j}=lasso(conf.Phi,conf.Phit,...
            %     conf.Psi,conf.Psit,conf.y,initSig,opt);
            % save(filename,'npg','-append');
            % 
            % u=10.^[-6 -6 -6 -6 -5 -5];
            % opt.u = u(i)*10^(j-2);
            % opt.continuation=false; opt.alphaStep='FISTA_L1';
            % npgs{i,j}=lasso(conf.Phi,conf.Phit,...
            %     conf.Psi,conf.Psit,conf.y,initSig,opt);
            % save(filename,'npgs','-append');

            u=10.^[-6 -6 -6 -6 -5 -4];
            opt.u = u(i)*10^(j-2);
            opt.noiseType='gaussian'; opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgLin{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,-log(conf.y),initSig,opt);
            save(filename,'npgLin','-append');
            opt.noiseType='poissonLogLink';
        end
    end
end

% the last example with glassbeads under the poisson model with log link 
% function
% Different from 009, the snr=1e4
if(any(runList==019))
    filename = [mfilename '_019.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    conf=ConfigCT();
    conf.imageName = 'glassBeadsSim';
    conf.PhiMode = 'gpuPrj';
    conf.PhiModeGen = 'gpuPrj';
    conf.dist = 17000;
    conf.beamharden = false;

    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    opt.maxItr=2e3; opt.thresh=1e-6; opt.snr=1e4; opt.debugLevel=1;
    opt.noiseType='poissonLogLink'; %'gaussian'; %
    for i=1:length(prjFull)-2
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull;
        opt=conf.setup(opt);
        conf.y = conf.y/max(conf.y(:)); % important to normalize
        initSig = maskFunc(conf.FBP(-log(conf.y)),opt.mask~=0);
        % initSig = opt.trueAlpha;

        for j=3:3
            fprintf('%s, i=%d, j=%d\n','X-ray CT example glassBeads Simulated',i,j);

            % fbp{i,j}.img=conf.FBP(-log(conf.y));
            % fbp{i,j}.alpha=fbp{i,j}.img(opt.mask~=0);
            % fbp{i,j}.RSE=sqrNorm(conf.y-conf.Phi(fbp{i,j}.alpha))/sqrNorm(conf.y);
            % fbp{i,j}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
            % fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);
            % save(filename,'fbp','-append');

            % u=10.^[-4 -4 -4 -4 -3 -3];
            % opt.u = u(i)*10^(j-2);
            % opt.continuation = false; opt.alphaStep='FISTA_ADMM_NNL1';
            % npg{i,j}=lasso(conf.Phi,conf.Phit,...
            %     conf.Psi,conf.Psit,conf.y,initSig,opt);
            % save(filename,'npg','-append');
             
            % u=10.^[-6 -6 -6 -6 -5 -5];
            % opt.u = u(i)*10^(j-2);
            % opt.continuation=false; opt.alphaStep='FISTA_L1';
            % npgs{i,j}=lasso(conf.Phi,conf.Phit,...
            %     conf.Psi,conf.Psit,conf.y,initSig,opt);
            % save(filename,'npgs','-append');

            u=10.^[-3 -3 -3 -3 -3 -3];
            opt.u = u(i)*10^(j-2);
            opt.noiseType='gaussian'; opt.continuation=false;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgLin{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,-log(conf.y),initSig,opt);
            save(filename,'npgLin','-append');
            opt.noiseType='poissonLogLink';
        end
    end
end

end
