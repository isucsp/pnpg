function [conf,opt] = runAsilomar2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0) runList = [0]; elseif(isempty(runList)) return; end
paperDir = '~/research/myPaper/asilomar2014/';

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
        npgsc{i}=Wrapper.NPGsc(conf.Phi,conf.Phit,conf.Psi,conf.Psit,y-gn.a0(i),zeros(p,1),opt);
        npgs{i}=Wrapper.NPGs(conf.Phi,conf.Phit,conf.Psi,conf.Psit,y-gn.a0(i),zeros(p,1),opt);
        fpcas{i}=Wrapper.FPCas(conf.Phi,conf.Phit,conf.Psi,conf.Psit,y-gn.a0(i),zeros(p,1),opt);
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
            fprintf('%s, i=%d, j=%d, k=%d\n','NPG',i,j,k);
            opt.u = u*10^(j-2);

            opt.continuation=false; opt.alphaStep='NPG';
            npg{i,j,k}=solver(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npg','-append');

            opt.continuation=false; opt.alphaStep='NPGs';
            npgs{i,j,k}=solver(conf.Phi,conf.Phit,...
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
            fprintf('%s, i=%d, j=%d\n','NPG',i,j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=false;
            opt.alphaStep='NPG';
            npg021{i,j}=solver(conf.Phi,conf.Phit,...
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
            fprintf('%s, i=%d, j=%d\n','NPG',i,j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=false;
            opt.alphaStep='NPG';
            npg031{i,j}=solver(conf.Phi,conf.Phit,...
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
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
    m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u = [1e-3,1e-3,1e-4,1e-4,1e-5,1e-5,1e-6,1e-6,1e-6];
    for k=1:5
        for i=1:length(m)
            opt.m=m(i); opt.snr=inf;
            [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
            initSig = Phit(invEAAt*y);

            opt.u = u(i)*10.^(-2:2);
            %gnet{i,k}=Wrapper.glmnet(Phi,wvltMat(length(opt.trueAlpha),dwt_L,daub),y,initSig,opt);

            for j=1:5
                fprintf('%s, i=%d, j=%d, k=%d\n','NPG',i,j,k);
                opt.u = u(i)*10^(j-3)*pNorm(Psit(Phit(y)),inf);

%               temp=opt; opt.thresh=1e-12; opt.maxItr=5e4;
%               % pgc12{i,j,k}=Wrapper.PGc(Phi,Phit,Psi,Psit,y,initSig,opt);
%               %sparsn12{i,j,k}=Wrapper.SpaRSAp(Phi,Phit,Psi,Psit,y,initSig,opt);
%               %spiral12{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
%               sparsa12 {i,j,k}=Wrapper.SpaRSA   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               opt=temp;
%               npgsc    {i,j,k}=Wrapper.NPGsc    (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgs     {i,j,k}=Wrapper.NPGs     (Phi,Phit,Psi,Psit,y,initSig,opt);

%               continue;

%               npgsT {i,j,k}=Wrapper.NPGs   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               temp=opt; opt.initStep='fixed';
%               fistal{i,j,k}=Wrapper.FISTA(Phi,Phit,Psi,Psit,y,initSig,opt);
%               opt=temp;
%               continue;

%               npgc     {i,j,k}=Wrapper.NPGc     (Phi,Phit,Psi,Psit,y,initSig,opt);
%               continue;

%               fpc      {i,j,k}=Wrapper.FPC      (Phi,Phit,Psi,Psit,y,initSig,opt);
%               sparsa   {i,j,k}=Wrapper.SpaRSA   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgc_nads{i,j,k}=Wrapper.NPGc_nads(Phi,Phit,Psi,Psit,y,initSig,opt);
%               npg_nads {i,j,k}=Wrapper.NPG_nads (Phi,Phit,Psi,Psit,y,initSig,opt);
%               pgc      {i,j,k}=Wrapper.PGc      (Phi,Phit,Psi,Psit,y,initSig,opt);
%               sparsn   {i,j,k}=Wrapper.SpaRSAp  (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npg      {i,j,k}=Wrapper.NPG     (Phi,Phit,Psi,Psit,y,initSig,opt);
%               spiral   {i,j,k}=Wrapper.SPIRAL   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               fista    {i,j,k}=Wrapper.FISTA    (Phi,Phit,Psi,Psit,y,initSig,opt);
                fpcas    {i,j,k}=Wrapper.FPCas    (Phi,Phit,Psi,Psit,y,initSig,opt);
            end

            save(filename);
        end
    end
end

% PET example
if(any(runList==003))
    filename = [mfilename '_003.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1; opt.noiseType='poisson';

    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    K=12;

    a  = [0, 0, 0, 0, 0, -1];
    as = [2, 2, 2, 1, 1,  1];
    aa = (3:-0.5:-6);

    for k=1:K
        for i=1:length(count)
            [y,Phi,Phit,Psi,Psit,fbpfunc,opt]=loadPET(count(i),opt);
            fbp{i,1,k}.alpha=maskFunc(fbpfunc(y),opt.mask~=0);
            fbp{i,1,k}.RMSE=sqrNorm(fbp{i,1,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);

            fprintf('fbp RMSE=%f\n',sqrNorm(fbp{i,1,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha));
            fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
            %u_max=pNorm(Psit(Phit(1-y./opt.bb(:))),inf);
            %u_max=pNorm(Psit(Phit(y)),inf);
            u_max=1;

            initSig=max(fbp{i,1,k}.alpha,0);
            %initSig=opt.trueAlpha;
            
            if(k>2) return; end

%           opt.fullcont=true;
%           opt.u=(10.^aa)*u_max; opt.maxItr=1e4; opt.thresh=1e-12;
%           npgFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgFull{i,k};
%           fprintf('k=%d, good a = 1e%g\n',k,max((aa(out.contRMSE==min(out.contRMSE)))));
%           opt.fullcont=false;

            j=1;
            fprintf('%s, i=%d, j=%d, k=%d\n','PET Example_003',i,j,k);
            opt.u = 10^a(i)*u_max;
%           npg   {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
%           spiral{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
%           npgc  {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
            if(i==5 && k==2)
                npgsc_s=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);
            end

            opt.u = 10^as(i)*u_max;
%           npgs  {i,j,k}=Wrapper.NPGs   (Phi,Phit,Psi,Psit,y,initSig,opt);
%           npgsc {i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);

%           % following are methods for weighted versions
%           ty=max(sqrt(y),1);
%           wPhi=@(xx) Phi(xx)./ty;
%           wPhit=@(xx) Phit(xx./ty);
%           wy=(y-opt.bb(:))./ty;
%           wu_max=pNorm(Psit(wPhit(wy)),inf);
%           opt.noiseType='gaussian';

%           opt.fullcont=true;
%           opt.u=(10.^aa)*wu_max; opt.maxItr=1e4; opt.thresh=1e-12;
%           wnpgFull {i,k}=Wrapper.NPG(wPhi,wPhit,Psi,Psit,wy,initSig,opt); out=wnpgFull{i,k};
%           fprintf('k=%d, good a = 1e%g\n',k,max((aa(out.contRMSE==min(out.contRMSE)))));
%           opt.fullcont=false;

%           opt.u = 10^a(i)*u_max;
%           fprintf('%s, i=%d, j=%d, k=%d\n','PET Example_003',i,1,k);
%           wnpg{i,k}=Wrapper.NPG         (wPhi,wPhit,Psi,Psit,wy,initSig,opt);
%           wspiral{i,k}=Wrapper.SPIRAL (wPhi,wPhit,Psi,Psit,wy,initSig,opt);
%           % wnpgc  {i,k}=Wrapper.NPGc   (wPhi,wPhit,Psi,Psit,wy,initSig,opt);
            save(filename);
        end
    end
end

% Linear model with Gaussian noise for Phantom example
% vary the number of measurements, and noise variance
if(any(runList==004))
    filename = [mfilename '_004.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    conf.imageName = 'phantom'; conf.imgSize=512;
    conf.beamharden = false;
    conf.PhiMode = 'gpuPrj';
    conf.PhiModeGen = 'gpuPrj';

    snr=[1e6 1e4 1e5 1e7];
    prjFull = [60, 80, 100, 120, 180, 360];
    a=[...
        -7, -7, -7, -7, -6, -6;...
        -5, -5, -5, -5, -5, -5;...
        -7, -7, -7, -6, -6, -6;...
        -7, -7, -7, -7, -7, -7;...
        ];
    aa = -1:-1:-8;  % old aa
    aa = -3:-1:-10;  % new aa

    opt.maxItr=1e4; opt.thresh=1e-6;
    opt.debugLevel=1;
    opt.noiseType='gaussian';

    K=2;
    for k=1:K
    for j=1:4
        opt.snr=snr(j);
        for i=1:length(prjFull)
            fprintf('%s, i=%d, j=%d, k=%d\n','Linear Gaussian Phantom Example: ',i,j,k);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
            opt=conf.setup(opt);
            initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
            initSig(initSig<0)=0;
            u_max=pNorm(conf.Psit(conf.Phit(conf.y)),inf);

            if(k>1) return; end

            if(j~=1 && (i~=3)) continue; end

%           temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
%           npgsFull{i,j}=Wrapper.NPGs(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           npgFull{i,j}=Wrapper.NPG(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           out=npgFull{i,j}; oRSE=[out.contRMSE(:); out.RMSE(end)];
%           fprintf('npg: i=%d, j=%d, good a = 1e%g\n',i,j,max(aa(oRSE==min(oRSE))));
%           out=npgsFull{i,j}; oRSE=[out.contRMSE(:); out.RMSE(end)];
%           fprintf('npgs: i=%d, j=%d, good a = 1e%g\n',i,j,max(aa(oRSE==min(oRSE))));
%           opt=temp;
%           save(filename); continue;

%           tic
%           fbp{i,j,k}.img=conf.FBP(conf.y);
%           fbp{i,j,k}.time=toc;
%           fbp{i,j,k}.alpha=fbp{i,j,k}.img(opt.mask~=0);
%           fbp{i,j,k}.RSE=sqrNorm(conf.y-conf.Phi(fbp{i,j,k}.alpha))/sqrNorm(conf.y);
%           fbp{i,j,k}.RMSE=sqrNorm(fbp{i,j,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
%           fprintf('fbp RMSE=%g\n',fbp{i,j,k}.RMSE);

            opt.u = 10^a(j,i)*u_max;
%           fpc   {i,j,k}=Wrapper.FPC    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           fista {i,j,k}=Wrapper.FISTA  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           fpcas {i,j,k}=Wrapper.FPCas  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

%           if(k>5)
%           fpcas {i,j,k}=Wrapper.FPCas  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           npg   {i,j,k}=Wrapper.NPG    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           end

%           npgs  {i,j,k}=Wrapper.NPGs   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           npgsSy{i,j,k}=Wrapper.NPGs_syn(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           npgc  {i,j,k}=Wrapper.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           npgsc {i,j,k}=Wrapper.NPGsc  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           sparsa{i,j,k}=Wrapper.SpaRSA (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           sparsn{i,j,k}=Wrapper.SpaRSAp(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           spiral{i,j,k}=Wrapper.SPIRAL (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

            if(k==1 && j==1 && i==3)
                temp=opt; opt.thresh=1e-12;

%               npgc_special1{i,j,k}=Wrapper.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
                npg12{i,j,k}=Wrapper.NPG   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%               spiral12{i,j,k}=Wrapper.SPIRAL   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%               sparsn12{i,j,k}=Wrapper.SpaRSAp  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%               sparsa12{i,j,k}=Wrapper.SpaRSA   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

                opt.u=opt.u*1e-4;

%               npgc_special2{i,j,k}=Wrapper.NPGc   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%               spiral_special2{i,j,k}=Wrapper.SPIRAL   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

                opt=temp;
            else
                continue;
            end

            save(filename);
        end
    end
    end
end

% vary the SNR for m=600, u is picked to give the best result
if(any(runList==005))
    filename = [mfilename '_005.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    opt.debugLevel=0; opt.maxItr=1e4; opt.thresh=1e-6;
    m  =[ 600 ];
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7];
    a  =[1e-1,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3,1e-4,1e-4,1e-5];
    aa =10.^(0:-1:-10);

    for k=1:5
        for i=1:length(snr)
            opt.m=m; opt.snr=snr(i);
            [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
            initSig = Phit(invEAAt*y);
            u_max=pNorm(Psit(Phit(y)),inf);

            % opt.fullcont=true;
            % opt.u=aa*u_max;
            % npgFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgFull{i,k};
            % fprintf('i=%d, good a = %e\n',i,aa(out.contRMSE==min(out.contRMSE)));
            % npgsFull{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgsFull{i,k};
            % fprintf('i=%d, good a = %e\n',i,aa(out.contRMSE==min(out.contRMSE)));
            % opt.fullcont=false;
            % save(filename); continue;

            j=1;
            opt.u = a(i)*u_max;
            fprintf('%s, i=%d, j=%d, k=%d\n','Example_005',i,j,k);
            npg   {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
            npgc  {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
            sparsn{i,j,k}=Wrapper.SpaRSAp(Phi,Phit,Psi,Psit,y,initSig,opt);
            spiral{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
            fista {i,j,k}=Wrapper.FISTA  (Phi,Phit,Psi,Psit,y,initSig,opt);
            npgsc {i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);
            fpcas {i,j,k}=Wrapper.FPCas  (Phi,Phit,Psi,Psit,y,initSig,opt);
            sparsa{i,j,k}=Wrapper.SpaRSA (Phi,Phit,Psi,Psit,y,initSig,opt);
            save(filename);
        end
    end
end

% Poisson example
% vary the number of measurements, with continuation, SNR=1e7
if(any(runList==006))
    filename = [mfilename '_006.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
    opt.noiseType='poisson'; opt.matrixType='nonneg'; opt.snr=1e7;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
    aa =10.^(3:-1:-10);

    for k=1:5
        for i=1:length(m)
            opt.m=m(i);
            [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
            initSig=randn(size(initSig));
            initSig=Phit(invEAAt*y)*0+1;
            fprintf('min=%d, max=%d\n',min(y), max(y));
            u_max=1;

            opt.trueAlpha=zeros(size(opt.trueAlpha)); y=y*0;
            opt.u=1e-3;

            Tnpg   =Wrapper.NPG   (Phi,Phit,Psi,Psit,y,initSig,opt);
            norm(Tnpg.alpha)

            Tspiral=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
            norm(Tspiral.alpha)

            keyboard

%           opt.fullcont=true;
%           opt.u=aa*u_max;
%           npgFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgFull{i,k};
%           fprintf('i=%d, good a = 1e%g\n',i,max(log10(aa(out.contRMSE==min(out.contRMSE)))));
%           npgsFull{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgsFull{i,k};
%           fprintf('i=%d, good a = 1e%g\n',i,max(log10(aa(out.contRMSE==min(out.contRMSE)))));
%           opt.fullcont=false;
%           save(filename); continue;

            for j=1:5;
                opt.u = a(i)*u_max*10^(j-3);

                fprintf('%s, i=%d, j=%d, k=%d\n','Example_006',i,j,k);
                npg    {i,j,k}=Wrapper.NPG   (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgc   {i,j,k}=Wrapper.NPGc  (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgs   {i,j,k}=Wrapper.NPGs  (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgsc  {i,j,k}=Wrapper.NPGsc (Phi,Phit,Psi,Psit,y,initSig,opt);
                spiral {i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
                temp=opt; opt.thresh=1e-8;
                spiral8{i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt=temp;
            end
            save(filename);
        end
    end
end

% Poisson example
% vary the number of measurements, with continuation, SNR=1e4
if(any(runList==016))
    filename = [mfilename '_016.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
    opt.noiseType='poisson2'; opt.matrixType='nonneg2'; opt.snr=1e7;
    m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400]; % should go from 200
    aa =10.^(3:-1:-10);

    for k=1:1
        for i=1:length(m)
            opt.m=m(i);
            [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
            initSig=Phit(invEAAt*y); initSig(initSig<0)=0;
            fprintf('min=%d, max=%d\n',min(y), max(y));
            u_max=1;

            opt.fullcont=true;
            opt.u=aa*u_max;
            npgFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgFull{i,k};
            fprintf('i=%d, good a = 1e%g\n',i,max(log10(aa(out.contRMSE==min(out.contRMSE)))));
            npgsFull{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgsFull{i,k};
            fprintf('i=%d, good a = 1e%g\n',i,max(log10(aa(out.contRMSE==min(out.contRMSE)))));
            opt.fullcont=false;
            save(filename); continue;

%           for j=1:5;
%               opt.u = a(i)*u_max*10^(j-2);
%               fprintf('%s, i=%d, j=%d, k=%d\n','Example_006',i,j,k);
%               npg   {i,j,k}=Wrapper.NPG   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgs  {i,j,k}=Wrapper.NPGs  (Phi,Phit,Psi,Psit,y,initSig,opt);
%               spiral{i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
%           end
            save(filename);
        end
    end
end

% skyline log link Poisson example
% vary the number of measurements, with continuation
if(any(runList==007))
    filename = [mfilename '_007.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    opt.maxItr=1e5; opt.thresh=1e-6;
    m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200
    aa=-1:-1:-8;
    aa=-3:-0.2:-6;
    opt.noiseType='poissonLogLink'; opt.matrixType='conv';
    for k=1:5
        for i=1:length(m)
            opt.m=m(i);
            [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
            initSig=-Phit(invEAAt*log(max(y,1)/max(y)))*0;
            fprintf('i=%d, k=%d, min=%d, max=%d\n',i,k,min(y), max(y));
            
            % if(any([2 3 4 6]==i)) continue; end

            %% fit by approximated Gaussian model for known I0
            % u_max=pNorm(Psit(Phit(y.*log(opt.I0./max(y,1)))),inf);
            %% fit by Poisson model with known I0
            % u_max=pNorm(Psit(Phit(y-opt.I0)),inf);
            %% fit by the unkown I0 log link Poisson model
            u_max=pNorm(Psit(Phit(y-mean(y))),inf);

            temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max; opt.noiseType='poissonLogLink';
%           gnet    {i,k}=Wrapper.glmnet(Utils.getMat(Phi,length(initSig)),Utils.getMat(Psi,length(Psit(initSig))),y,initSig,opt);
%           npgFull {i,k}=Wrapper.NPG   (Phi,Phit,Psi,Psit,y,initSig,opt);
%           npgsFull{i,k}=Wrapper.NPGs  (Phi,Phit,Psi,Psit,y,initSig,opt);
            initSig=-Phit(invEAAt*log(max(y,1)/max(y)));
            npgsFull1{i,k}=Wrapper.NPGs  (Phi,Phit,Psi,Psit,y,initSig,opt);
            opt=temp;
            continue;

            temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max; opt.noiseType='poissonLogLink0';
            npgFull_knownI0 {i,k}=Wrapper.NPG (Phi,Phit,Psi,Psit,y,initSig,opt);
            npgsFull_knownI0{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt);
            opt=temp;
           
            %% fit by approximated Gaussian model for known I0
            yy=log(opt.I0./max(y,1)); yy=yy.*sqrt(y);
            wPhi=@(xxx) sqrt(y).*Phi(xxx); wPhit=@(xxx) Phit(sqrt(y).*xxx);

            temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max; opt.noiseType='gaussian';
            npglwFull {i,k}=Wrapper.NPG(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
            npgslwFull{i,k}=Wrapper.NPGs(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
            opt=temp;

            save(filename);
            continue

            yy=log(opt.I0./max(y,1));
            u_max=pNorm(Psit(Phit(yy)),inf);
            temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max; opt.noiseType='gaussian';
            npglFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,yy,initSig,opt);
            npgslFull{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,yy,initSig,opt);
            opt=temp;
            continue;

            %% fit by approximated Gaussian model for unknown I0

            opt.u = u(i)*10^(j-3)*u_max;
            fprintf('%s, i=%d, j=%d, k=%d\n','NPG',i,j,k);

            opt.noiseType='gaussian';
            temp=y; temp(temp==0)=1; temp=log(opt.I0./temp); temp=temp.*sqrt(y);
            wPhi=@(xxx) sqrt(y).*Phi(xxx);
            wPhit=@(xxx) Phit(sqrt(y).*xxx);
            npgclw {i,j,k}=Wrapper.NPGc   (wPhi,wPhit,Psi,Psit,temp,initSig,opt);
            npgsclw{i,j,k}=Wrapper.NPGsc  (wPhi,wPhit,Psi,Psit,temp,initSig,opt);

            temp=y; temp(temp==0)=1; temp=log(opt.I0./temp);
            npgcl {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,temp,initSig,opt);
            npgscl{i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,temp,initSig,opt);
            opt.noiseType='poissonLogLink';

            npgc  {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
            npgsc {i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);
        end
        save(filename);
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

        sparsn   {i,j}=Wrapper.SpaRSAp  (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename);
        continue
        npgc_nads{i,j}=Wrapper.NPGc_nads(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npg_nads {i,j}=Wrapper.NPG_nads (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npg      {i,j}=Wrapper.NPG     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npgs     {i,j}=Wrapper.NPGs     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

        npgsc    {i,j}=Wrapper.NPGsc    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        npgc     {i,j}=Wrapper.NPGc     (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);

        spiral   {i,j}=Wrapper.SPIRAL   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        fpcas    {i,j}=Wrapper.FPCas    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        sparsa   {i,j}=Wrapper.SpaRSA   (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
        fista    {i,j}=Wrapper.FISTA    (conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
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
    conf.PhiMode = 'cpuPrj';    % change the option to cpuPrj if no GPU equipped
    conf.PhiModeGen = 'cpuPrj'; % change the option to cpuPrj if no GPU equipped
    conf.dist = 17000;
    conf.beamharden = false;

    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    aa = -1:-1:-10;
    opt.maxItr=4e3; opt.thresh=1e-6; opt.snr=1e6; opt.debugLevel=1;
    opt.noiseType='poissonLogLink'; %'gaussian'; %
    for i=[5 6 1 2 3 4]
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0)); j=1;
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull;
        opt=conf.setup(opt);
        fprintf('i=%d, min=%d, max=%d\n',i,min(conf.y), max(conf.y));
        initSig = maskFunc(conf.FBP(-log(max(conf.y,1)/max(conf.y))),opt.mask~=0);

        fbp{i,j}.img=conf.FBP(-log(max(conf.y,1)/max(conf.y)));
        fbp{i,j}.alpha=fbp{i,j}.img(opt.mask~=0);
        fbp{i,j}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);
        fprintf('fbp after truncation RMSE=%g\n',rmseTruncate(fbp{i,j},opt.trueAlpha));

        if(i~=2) continue; end
        % the poisson model with log link, where I0 is unknown
        % initSig = opt.trueAlpha;
        % u_max=pNorm(conf.Psit(conf.Phit(conf.y-opt.I0)),inf); % for loglink0
        % u_max=pNorm(conf.Psit(conf.Phit(conf.y.*log(max(conf.y,1)/opt.I0))),inf) % for loglink0 approximated by weight
        u_max=pNorm(conf.Psit(conf.Phit(conf.y-mean(conf.y))),inf); % for loglink

%       opt.fullcont=true;
%       opt.u=10.^aa*u_max;
%       npgFull{i}=Wrapper.NPG(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt); out=npgFull{i};
%       fprintf('i=%d, good a = 1e%g\n',i,max((aa(out.contRMSE==min(out.contRMSE)))));
%       npgsFull{i}=Wrapper.NPGs(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt); out=npgsFull{i};
%       fprintf('i=%d, good a = 1e%g\n',i,max((aa(out.contRMSE==min(out.contRMSE)))));
%       opt.fullcont=false;
%       save(filename); continue;

%       a=[-3.5 -3.75 -4 -4.25 -4.5];
%       for j=1:5;
%           fprintf('%s, i=%d, j=%d\n','X-ray CT example glassBeads Simulated',i,j);
%           opt.u = 10^a(j)*u_max;
%           npg{i,j}=Wrapper.NPGc(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           npgs{i,j}=Wrapper.NPGsc(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt);
%           save(filename);
%       end
%       continue

        % fit with the poisson model with log link but known I0
        temp=opt; opt.I0=conf.I0;

        opt.noiseType='poissonLogLink0';
        opt.fullcont=true;
        opt.u=10.^aa*u_max;
        npg0Full{i}=Wrapper.NPG(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt); out=npg0Full{i};
        fprintf('i=%d, good a = 1e%g\n',i,max((aa(out.contRMSE==min(out.contRMSE)))));
        npgs0Full{i}=Wrapper.NPGs(conf.Phi,conf.Phit,conf.Psi,conf.Psit,conf.y,initSig,opt); out=npgs0Full{i};
        fprintf('i=%d, good a = 1e%g\n',i,max((aa(out.contRMSE==min(out.contRMSE)))));
        opt.fullcont=false;

%       opt.noiseType='gaussian';
%       wPhi=@(xx) sqrt(conf.y).*conf.Phi(xx);
%       wPhit=@(xx) conf.Phit(sqrt(conf.y).*xx);
%       wy=sqrt(conf.y).*(log(opt.I0)-log(max(conf.y,1)));
%       u_max=pNorm(conf.Psit(wPhit(wy)),inf);
%       opt.fullcont=true;
%       opt.u=10.^aa*u_max;
%       wnpgFull{i}=Wrapper.NPG(wPhi,wPhit,conf.Psi,conf.Psit,wy,initSig,opt); out=wnpgFull{i};
%       fprintf('i=%d, good a = 1e%g\n',i,max((aa(out.contRMSE==min(out.contRMSE)))));
%       wnpgsFull{i}=Wrapper.NPGs(wPhi,wPhit,conf.Psi,conf.Psit,wy,initSig,opt); out=wnpgsFull{i};
%       fprintf('i=%d, good a = 1e%g\n',i,max((aa(out.contRMSE==min(out.contRMSE)))));
%       opt.fullcont=false;

        opt=temp;
         
        % fit with the approximated Gaussian model

        save(filename);
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
            % opt.continuation = false; opt.alphaStep='NPG';
            % npg{i,j}=solver(conf.Phi,conf.Phit,...
            %     conf.Psi,conf.Psit,conf.y,initSig,opt);
            % save(filename,'npg','-append');
             
            % u=10.^[-6 -6 -6 -6 -5 -5];
            % opt.u = u(i)*10^(j-2);
            % opt.continuation=false; opt.alphaStep='NPGs';
            % npgs{i,j}=solver(conf.Phi,conf.Phit,...
            %     conf.Psi,conf.Psit,conf.y,initSig,opt);
            % save(filename,'npgs','-append');

            u=10.^[-3 -3 -3 -3 -3 -3];
            opt.u = u(i)*10^(j-2);
            opt.noiseType='gaussian'; opt.continuation=false;
            opt.alphaStep='NPG';
            npgLin{i,j}=solver(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,-log(conf.y),initSig,opt);
            save(filename,'npgLin','-append');
            opt.noiseType='poissonLogLink';
        end
    end
end

if(any(runList==902))
    filename = [mfilename '_002.mat'];
    load(filename);

    m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u = [1e-3,1e-3,1e-4,1e-4,1e-5,1e-5,1e-6,1e-6,1e-6];
    idx=2:2:7;
    K = 1;

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    npgscTime = mean(Cell.getField( npgsc(:,:,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,:,1:K),'time'),3);
    fpcasTime = mean(Cell.getField( fpcas(:,:,1:K),'cpu' ),3);
    fpcTime   = mean(Cell.getField(   fpc(:,:,1:K),'time' ),3);
    fistaTime = mean(Cell.getField( fista(:,:,1:K),'time'),3);
    sparsaTime= mean(Cell.getField(sparsa(:,:,1:K),'time'),3);
    sparsnTime= mean(Cell.getField(sparsn(:,:,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    npgscCost = mean(Cell.getField( npgsc(:,:,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,:,1:K),'cost'),3);
    fpcasCost = mean(Cell.getField( fpcas(:,:,1:K),'f'   ),3);
    fpcCost   = mean(Cell.getField(   fpc(:,:,1:K),'cost' ),3);
    fistaCost = mean(Cell.getField( fista(:,:,1:K),'cost'),3);
    sparsaCost= mean(Cell.getField(sparsa(:,:,1:K),'cost'),3);
    sparsnCost= mean(Cell.getField(sparsn(:,:,1:K),'cost'),3);

    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,:,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,:,1:K),'reconerror'),3);
    fpcasRMSE = mean(Cell.getField( fpcas(:,:,1:K),'RMSE'),3);
    fpcRMSE   = mean(Cell.getField(   fpc(:,:,1:K),'RMSE' ),3);
    fistaRMSE = mean(Cell.getField( fista(:,:,1:K),'RMSE'),3);
    sparsaRMSE= mean(Cell.getField(sparsa(:,:,1:K),'RMSE'),3);
    sparsnRMSE= mean(Cell.getField(sparsn(:,:,1:K),'RMSE'),3);

    npgnasTime = mean(Cell.getField(   npg_nads(:,:,1:K),'time'),3);
    npgcnasTime= mean(Cell.getField(  npgc_nads(:,:,1:K),'time'),3);
    npgnasCost = mean(Cell.getField(   npg_nads(:,:,1:K),'cost'),3);
    npgcnasCost= mean(Cell.getField(  npgc_nads(:,:,1:K),'cost'),3);
    npgnasRMSE = mean(Cell.getField(   npg_nads(:,:,1:K),'RMSE'),3);
    npgcnasRMSE= mean(Cell.getField(  npgc_nads(:,:,1:K),'RMSE'),3);

    for i=1:length(m)
        temp=[];
        for k=1:K
            temp(:,k)=gnet{i,k}.RMSE(:);
        end
%       gnetRMSE(i,:)=mean(temp,2)';
    end

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c4]=find(spiralRMSE== repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c8]=find(sparsnRMSE== repmat(min(sparsnRMSE,[],2),1,5)); [r,idx8]=sort(r);

    [r,c3]=find(  npgsRMSE== repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r);
    [r,c7]=find(sparsaRMSE== repmat(min(sparsaRMSE,[],2),1,5)); [r,idx7]=sort(r);
    [r,c9]=find( npgscRMSE== repmat(min( npgscRMSE,[],2),1,5)); [r,idx9]=sort(r);
    disp([c1(idx1), c2(idx2), c4(idx4), c8(idx8) zeros(9,1) c3(idx3), c5(idx5), c6(idx6), c7(idx7), c9(idx9) ]);
    keyboard
    uNonneg=[3 3 3 3 4 4 4 4 3];
       uNeg=[4 4 4 4 4 4 4 4 3];
    figure;
    semilogy(m,   npgRMSE((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcRMSE((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m, fpcasRMSE((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaRMSE((c6(idx6)-1)*9+(1:9)'),'b-.');
    semilogy(m,sparsaRMSE((c7(idx7)-1)*9+(1:9)'),'y-p');
    semilogy(m,sparsnRMSE((c8(idx8)-1)*9+(1:9)'),'r-x');
    semilogy(m, npgscRMSE((c9(idx9)-1)*9+(1:9)'),'k-.');
    legend('npg','npgc','npgs','spiral','fpcas','fista','sparsa','sparsa','npgsc');
    figure;
    semilogy(m,   npgTime((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcTime((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m, fpcasTime((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaTime((c6(idx6)-1)*9+(1:9)'),'b-.');
    semilogy(m,sparsaTime((c7(idx7)-1)*9+(1:9)'),'y-p');
    semilogy(m,sparsnTime((c8(idx8)-1)*9+(1:9)'),'r-x');
    semilogy(m, npgscTime((c9(idx9)-1)*9+(1:9)'),'k-.');
    legend('npg','npgc','npgs','spiral','fpcas','fista','sparsa','sparsa','npgsc');
    figure;
    semilogy(m,   npgRMSE((uNonneg-1)*9+(1:9)),'r-*'); hold on;
    semilogy(m,  npgcRMSE((uNonneg-1)*9+(1:9)),'c-p');
    semilogy(m,spiralRMSE((uNonneg-1)*9+(1:9)),'k-^');
    semilogy(m,sparsnRMSE((uNonneg-1)*9+(1:9)),'r-x');
    semilogy(m,  npgsRMSE((uNonneg-1)*9+(1:9)),'k-s');
    semilogy(m, fpcasRMSE((uNonneg-1)*9+(1:9)),'g-o');
    semilogy(m, fistaRMSE((uNonneg-1)*9+(1:9)),'b-.');
    semilogy(m,sparsaRMSE((uNonneg-1)*9+(1:9)),'y-p');
%   semilogy(m,  gnetRMSE((uNonneg-1)*9+(1:9)),'r:>');
    semilogy(m, npgscRMSE((uNonneg-1)*9+(1:9)),'k-.');
    legend('npg','npgc','spiral','sparsa','npgs','fpcas','fista','sparsa','npgsc');
    figure;
    semilogy(m,   npgTime((uNonneg-1)*9+(1:9)),'r-*'); hold on;
    semilogy(m,  npgcTime((uNonneg-1)*9+(1:9)),'c-p');
    semilogy(m,spiralTime((uNonneg-1)*9+(1:9)),'k-^');
    semilogy(m,sparsnTime((uNonneg-1)*9+(1:9)),'r-x');
    semilogy(m,  npgsTime((uNonneg-1)*9+(1:9)),'k-s');
    semilogy(m, fpcasTime((uNonneg-1)*9+(1:9)),'g-o');
    semilogy(m, fistaTime((uNonneg-1)*9+(1:9)),'b-.');
    semilogy(m,sparsaTime((uNonneg-1)*9+(1:9)),'y-p');
    semilogy(m, npgscTime((uNonneg-1)*9+(1:9)),'k-.');
    legend('npg','npgc','spiral','sparsa','npgs','fpcas','fista','sparsa','npgsc');

    f=fopen('selectedTime.data','w');
    for mIdx=1:length(m)
        fprintf(f,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%s\t%s\n',...
                npgTime(mIdx,uNonneg(mIdx)), ...
               npgcTime(mIdx,uNonneg(mIdx)), ...
             spiralTime(mIdx,uNonneg(mIdx)), ...
             sparsnTime(mIdx,uNonneg(mIdx)), ...
               npgsTime(mIdx,uNonneg(mIdx)), ...
              npgscTime(mIdx,uNonneg(mIdx)), ...
              fpcasTime(mIdx,uNonneg(mIdx)), ...
              fistaTime(mIdx,uNonneg(mIdx)), ...
             sparsaTime(mIdx,uNonneg(mIdx)), ...
            m(mIdx),num2str(log10(u((mIdx)))+uNonneg(mIdx)-3), num2str(log10(u(mIdx))+uNeg(mIdx)-3));
    end
    fclose(f);

    keyboard

    as=1:5;
    forSave=[]; forTime=[];
    for mIdx=idx
        figure(900);
        semilogy(log10(u(mIdx))+as-3,    npgRMSE(mIdx,as),'r-*'); hold on;
        semilogy(log10(u(mIdx))+as-3,   npgcRMSE(mIdx,as),'r.-');
        semilogy(log10(u(mIdx))+as-3, sparsnRMSE(mIdx,as),'r-s');
        semilogy(log10(u(mIdx))+as-3, spiralRMSE(mIdx,as),'r-^');
        semilogy(log10(u(mIdx))+as-3,   npgsRMSE(mIdx,as),'k-s');
        semilogy(log10(u(mIdx))+as-3,  fpcasRMSE(mIdx,as),'g-o');
        semilogy(log10(u(mIdx))+as-3,  fistaRMSE(mIdx,as),'g-.');
        semilogy(log10(u(mIdx))+as-3, sparsaRMSE(mIdx,as),'g->');
        semilogy(log10(u(mIdx))+as-3,  npgscRMSE(mIdx,as),'g-*');
        semilogy(log10(u(mIdx))+as-3,    fpcRMSE(mIdx,as),'g:p');
        semilogy(gnet{mIdx,1}.a,   gnet{mIdx,1}.RMSE(:),'r:>');

        forSave=[forSave log10(u(mIdx))+as(:)-3];
        forSave=[forSave reshape(   npgRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(  npgcRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(sparsnRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(spiralRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(  npgsRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape( fpcasRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape( fistaRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(sparsaRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape( npgscRMSE(mIdx,as),[],1)];

        figure;
        semilogy(log10(u(mIdx))+as-3,    npgTime(mIdx,as),'r-*'); hold on;
        semilogy(log10(u(mIdx))+as-3,   npgcTime(mIdx,as),'r.-');
        semilogy(log10(u(mIdx))+as-3, sparsnTime(mIdx,as),'r-s');
        semilogy(log10(u(mIdx))+as-3, spiralTime(mIdx,as),'r-^');
        semilogy(log10(u(mIdx))+as-3,   npgsTime(mIdx,as),'k-s');
        semilogy(log10(u(mIdx))+as-3,  fpcasTime(mIdx,as),'g-o');
        semilogy(log10(u(mIdx))+as-3,  fistaTime(mIdx,as),'g-.');
        semilogy(log10(u(mIdx))+as-3, sparsaTime(mIdx,as),'g->');
        semilogy(log10(u(mIdx))+as-3,  npgscTime(mIdx,as),'g-*');
        semilogy(log10(u(mIdx))+as-3,    fpcTime(mIdx,as),'g:p');
        legend('npg','npgc','sparsn','spiral','npgs','fpcas','fista','sparsa','npgsc','fpc');
        title(sprintf('mIdx=%d',mIdx));

        forTime=[forTime log10(u(mIdx))+as(:)-3];
        forTime=[forTime reshape(   npgTime(mIdx,as),[],1)];
        forTime=[forTime reshape(  npgcTime(mIdx,as),[],1)];
        forTime=[forTime reshape(sparsnTime(mIdx,as),[],1)];
        forTime=[forTime reshape(spiralTime(mIdx,as),[],1)];
        forTime=[forTime reshape(  npgsTime(mIdx,as),[],1)];
        forTime=[forTime reshape( fpcasTime(mIdx,as),[],1)];
        forTime=[forTime reshape( fistaTime(mIdx,as),[],1)];
        forTime=[forTime reshape(sparsaTime(mIdx,as),[],1)];
        forTime=[forTime reshape( npgscTime(mIdx,as),[],1)];
    end
    figure(900); 
    legend('npg','npgc','sparsn','spiral','npgs','fpcas','fista','sparsa','npgsc','fpc','glmnet');
    save('rmseVsA.data','forSave','-ascii');
    save('timeVsA.data','forTime','-ascii');

    keyboard

    mIdx=6; as=gEle(c2(idx2),mIdx); forSave=[]; t=0;
    q=(1:max(find(sparsn12{mIdx,as}.cost(:)>=sparsn{mIdx,as}.cost(end))))';
    t=t+1; temp=      npg{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp; % 3rd col
    t=t+1; temp=      npg{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp; % 7th col
    t=t+1; temp=     npgc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp; % 11th col
    t=t+1; temp= sparsn12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    q=(1:max(find(spiral12{mIdx,as}.cost(:)>=spiral{mIdx,as}.cost(end))))';
    t=t+1; temp= spiral12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp; % 17th col
    t=t+1; temp= spiral12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    q=(1:max(find(sparsa12{mIdx,as}.cost(:)>=sparsa{mIdx,as}.cost(end))))';
    t=t+1; temp= sparsa12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=    fista{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp; % 41st col
    t=t+1; temp=    fista{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    fista{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    fista{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp; % 45th col
    t=t+1; temp= spiral12{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    save('traceLinGauss.data','forSave','-ascii');

    mc=forSave(:,[3,7,11,19]); mc = min(mc(mc(:)>0));
    figure;
    loglog(forSave(:, 2),forSave(:, 3)-mc,'r-'); hold on;
    loglog(forSave(:, 6),forSave(:, 7)-mc,'r-.');
    loglog(forSave(:,14),forSave(:,15)-mc,'c--');
    loglog(forSave(:,18),forSave(:,19)-mc,'b:'); 
    legend('npg','npgc','sprsa 12','spiral');

    mc=forSave(:,[23,27,31,35,39,43]); mc = min(mc(mc(:)>0));
    figure;
    loglog(forSave(:,22),forSave(:,23)-mc,'r-'); hold on;
    loglog(forSave(:,26),forSave(:,27)-mc,'g-.');
    loglog(forSave(:,30),forSave(:,31)-mc,'c--');
    loglog(forSave(:,34),forSave(:,35)-mc,'b:'); 
    loglog(forSave(:,38),forSave(:,39)-mc,'k-'); 
    loglog(forSave(:,42),forSave(:,43)-mc,'r--'); 
    legend('npgs','npgsc','sparsa','fpc','sparsa12','fista');

    keyboard

    mIdx=6; as=gEle(c2(idx2),mIdx); forSave=[]; t=0;
    t=t+1; temp=  npgc{mIdx,as}.RMSE(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.time(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.cost(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.difAlpha(:);  forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.uRecord(:,2); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.contThresh(:);forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.RMSE(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.time(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.cost(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.difAlpha(:);  forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.uRecord(:,2); forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.contThresh(:);forSave(1:length(temp),t)=temp;
    save('continuation.data','forSave','-ascii');

    keyboard

    temp = 4;
    signal=npg{1}.opt.trueAlpha;
    signal=[signal,    npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,   npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,   npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  fista{gEle((c6(idx6)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, sparsa{gEle((c7(idx7)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, sparsn{gEle((c8(idx8)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  npgsc{gEle((c9(idx9)-1)*9+(1:9)',temp)}.alpha];
    save('skyline.data','signal','-ascii');

[gEle((c1(idx1)-1),temp);
gEle((c2(idx2)-1),temp);
gEle((c3(idx3)-1),temp);
gEle((c4(idx4)-1),temp);
gEle((c5(idx5)-1),temp);
gEle((c6(idx6)-1),temp);
gEle((c7(idx7)-1),temp);
gEle((c8(idx8)-1),temp);
gEle((c9(idx9)-1),temp);]'

    figure; plot(signal(:,2)); hold on; plot(signal(:,1),'r'); title('NPG');
    figure; plot(signal(:,4)); hold on; plot(signal(:,1),'r'); title('NPGs');
    figure; plot(signal(:,6)); hold on; plot(signal(:,1),'r'); title('FPCas');
    fprintf('\nfor N=350:\n'); temp=4;
    fprintf('   npgRec RMSE: %g%% -> %g%%\n',   npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(   npg{gEle((c1(idx1)-1)*9+(1:9)',temp)})*100);
    fprintf('  npgcRec RMSE: %g%% -> %g%%\n',  npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(  npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)})*100);
    fprintf('  npgsRec RMSE: %g%% -> %g%%\n',  npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(  npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)})*100);
    fprintf('spiralRec RMSE: %g%% -> %g%%\n',spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)}.reconerror(end)*100,rmseTruncate(spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)})*100);
    fprintf(' fpcasRec RMSE: %g%% -> %g%%\n', fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate( fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)})*100);
    fprintf(' fistaRec RMSE: %g%% -> %g%%\n', fista{gEle((c6(idx6)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate( fista{gEle((c6(idx6)-1)*9+(1:9)',temp)})*100);

    fprintf('\nfor N=250:\n'); temp=2;
    fprintf('   npgRec RMSE: %g%% -> %g%%\n',   npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(   npg{gEle((c1(idx1)-1)*9+(1:9)',temp)})*100);
    fprintf('  npgcRec RMSE: %g%% -> %g%%\n',  npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(  npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)})*100);
    fprintf('  npgsRec RMSE: %g%% -> %g%%\n',  npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(  npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)})*100);
    fprintf('spiralRec RMSE: %g%% -> %g%%\n',spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)}.reconerror(end)*100,rmseTruncate(spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)})*100);
    fprintf(' fpcasRec RMSE: %g%% -> %g%%\n', fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate( fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)})*100);
    fprintf(' fistaRec RMSE: %g%% -> %g%%\n', fista{gEle((c6(idx6)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate( fista{gEle((c6(idx6)-1)*9+(1:9)',temp)})*100);

    fprintf('\nfor N=500:\n'); temp=6;
    fprintf('   npgRec RMSE: %g%% -> %g%%\n',   npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(   npg{gEle((c1(idx1)-1)*9+(1:9)',temp)})*100);
    fprintf('  npgcRec RMSE: %g%% -> %g%%\n',  npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(  npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)})*100);
    fprintf('  npgsRec RMSE: %g%% -> %g%%\n',  npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate(  npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)})*100);
    fprintf('spiralRec RMSE: %g%% -> %g%%\n',spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)}.reconerror(end)*100,rmseTruncate(spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)})*100);
    fprintf(' fpcasRec RMSE: %g%% -> %g%%\n', fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate( fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)})*100);
    fprintf(' fistaRec RMSE: %g%% -> %g%%\n', fista{gEle((c6(idx6)-1)*9+(1:9)',temp)}.RMSE(end)*100      ,rmseTruncate( fista{gEle((c6(idx6)-1)*9+(1:9)',temp)})*100);

    M=length(m);
    str=        '$m$            ';              for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%10d',str,m(i)); end; end;
    str=sprintf('%s\\\\\\hline',str);
    str=sprintf('%s\\\\\nNPG            ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,   npg{gEle((c1(idx1)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nNPG$_\\text{C}$ ',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,  npgc{gEle((c2(idx2)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nNPG$_\\text{S}$ ',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,  npgs{gEle((c3(idx3)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nSPIRAL         ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,spiral{gEle((c4(idx4)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nFPC$_\\text{AS}$',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str, fpcas{gEle((c5(idx5)-1)*9+(1:9)',i)}.f   (end));end; end;
    str=sprintf('%s\nFISTA          ', str);    for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str, fista{gEle((c6(idx6)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nSpaRSA         ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,spiral{gEle((c7(idx7)-1)*9+(1:9)',i)}.cost(end));end; end;
    file=fopen('varyMeasurementTable.tex','w'); fprintf(file,'%s',str); fclose(file);

    % figure;
    % for i=1:M;
    %     semilogy(npgs{i,idx3,1}.stepSize); hold on; semilogy(fista{i,idx6,1}.stepSize,'r:');
    %     semilogy([1,length(fista{i,idx6,1}.RMSE)],ones(1,2)*1/fista{i,idx6,1}.opt.L,'k-.');
    %     hold off;
    %     pause;
    % end
      npgItr=[];   
     npgcItr=[];
     npgsItr=[];
   spiralItr=[];
    fpcasItr=[];
    fistaItr=[];
   sparsaItr=[];
   sparsnItr=[];
    npgscItr=[];

    for i=1:K
        temp=   npg(:,:,i); temp=temp((c1(idx1)-1)*9+(1:9)');    npgItr=[   npgItr,showResult(temp,2,'p'   )];
        temp=  npgc(:,:,i); temp=temp((c2(idx2)-1)*9+(1:9)');   npgcItr=[  npgcItr,showResult(temp,2,'p'   )];
        temp=  npgs(:,:,i); temp=temp((c3(idx3)-1)*9+(1:9)');   npgsItr=[  npgsItr,showResult(temp,2,'p'   )];
        temp=spiral(:,:,i); temp=temp((c4(idx4)-1)*9+(1:9)'); spiralItr=[spiralItr,showResult(temp,2,'p'   )];
        temp= fpcas(:,:,i); temp=temp((c5(idx5)-1)*9+(1:9)');  fpcasItr=[ fpcasItr,showResult(temp,2,'itr' )];
        temp= fista(:,:,i); temp=temp((c6(idx6)-1)*9+(1:9)');  fistaItr=[ fistaItr,showResult(temp,2,'p'   )];
        temp=sparsa(:,:,i); temp=temp((c7(idx7)-1)*9+(1:9)'); sparsaItr=[sparsaItr,showResult(temp,3,'RMSE')];
        temp=sparsn(:,:,i); temp=temp((c8(idx8)-1)*9+(1:9)'); sparsnItr=[sparsnItr,showResult(temp,3,'RMSE')];
        temp= npgsc(:,:,i); temp=temp((c9(idx9)-1)*9+(1:9)');  npgscItr=[ npgscItr,showResult(temp,2,'p'   )];
    end

    forSave=[];
    forSave=[forSave,    npgTime((c1(idx1)-1)*9+(1:9)')];
    forSave=[forSave,   npgcTime((c2(idx2)-1)*9+(1:9)')];
    forSave=[forSave,   npgsTime((c3(idx3)-1)*9+(1:9)')];
    forSave=[forSave, spiralTime((c4(idx4)-1)*9+(1:9)')];
    forSave=[forSave,  fpcasTime((c5(idx5)-1)*9+(1:9)')];
    forSave=[forSave,  fistaTime((c6(idx6)-1)*9+(1:9)')];

    forSave=[forSave,    npgCost((c1(idx1)-1)*9+(1:9)')];
    forSave=[forSave,   npgcCost((c2(idx2)-1)*9+(1:9)')];
    forSave=[forSave,   npgsCost((c3(idx3)-1)*9+(1:9)')];
    forSave=[forSave, spiralCost((c4(idx4)-1)*9+(1:9)')];
    forSave=[forSave,  fpcasCost((c5(idx5)-1)*9+(1:9)')];
    forSave=[forSave,  fistaCost((c6(idx6)-1)*9+(1:9)')];

    forSave=[forSave,    npgRMSE((c1(idx1)-1)*9+(1:9)')];
    forSave=[forSave,   npgcRMSE((c2(idx2)-1)*9+(1:9)')];
    forSave=[forSave,   npgsRMSE((c3(idx3)-1)*9+(1:9)')];
    forSave=[forSave, spiralRMSE((c4(idx4)-1)*9+(1:9)')];
    forSave=[forSave,  fpcasRMSE((c5(idx5)-1)*9+(1:9)')];
    forSave=[forSave,  fistaRMSE((c6(idx6)-1)*9+(1:9)')];
    forSave=[forSave, m(:)];
    forSave=[forSave, sparsaTime((c7(idx7)-1)*9+(1:9)')];
    forSave=[forSave, sparsaCost((c7(idx7)-1)*9+(1:9)')];
    forSave=[forSave, sparsaRMSE((c7(idx7)-1)*9+(1:9)')];
    forSave=[forSave, sparsnTime((c8(idx8)-1)*9+(1:9)')];
    forSave=[forSave, sparsnCost((c8(idx8)-1)*9+(1:9)')];
    forSave=[forSave, sparsnRMSE((c8(idx8)-1)*9+(1:9)')];
    forSave=[forSave,  npgscTime((c9(idx9)-1)*9+(1:9)')];
    forSave=[forSave,  npgscCost((c9(idx9)-1)*9+(1:9)')];
    forSave=[forSave,  npgscRMSE((c9(idx9)-1)*9+(1:9)')];
    save('varyMeasurement.data','forSave','-ascii');

    forSave=m(:);
    forSave=[forSave,    npgTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,   npgcTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave, spiralTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave, sparsnTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,   npgsTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,  fpcasTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,  fistaTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave, sparsaTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,  npgscTime((uNonneg-1)*9+(1:9))'];
    save('varyMeasurementTime.data','forSave','-ascii');

    keyboard

    mIdx=2; experi=1; forSave=[]; t=0;
    npgsT=npgsT(:,:,experi); npgsn20T=npgs(:,:,experi); fistaT=fista(:,:,experi); fistalT=fistal(:,:,experi); fistalT{9,6}=[];
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c3(idx3),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c3(idx3),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c3(idx3),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c3(idx3),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c3(idx3),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c3(idx3),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c3(idx3),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c3(idx3),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    disp([ c3(idx3) c6(idx6)]);
    disp([   npgsT{mIdx,gEle(c3(idx3),mIdx)}.p;...
           fistalT{mIdx,gEle(c3(idx3),mIdx)}.p;...
            fistaT{mIdx,gEle(c6(idx6),mIdx)}.p]);
    disp([   npgsT{mIdx,gEle(c3(idx3),mIdx)}.time(end); ...
           fistalT{mIdx,gEle(c3(idx3),mIdx)}.time(end); ...
            fistaT{mIdx,gEle(c6(idx6),mIdx)}.time(end)]);
    temp=forSave(:,13:16); temp=temp(:); temp=temp(temp>0); temp=min(temp); forSave(:,13:16)=forSave(:,13:16)-temp;
    save('stepSizeLin.data','forSave','-ascii');
    figure(1); hold off; semilogy(forSave(:,9),forSave(:,5),'r'); hold on;
    semilogy(forSave(:,10),forSave(:,6),'g');
    semilogy(forSave(:,11),forSave(:,7),'b');
    semilogy(forSave(:,12),forSave(:,8),'c');
    legend('npgs','npgs20','fistaBB','fistaL');
    figure(2); hold off; semilogy(forSave(:,9),forSave(:,13),'r'); hold on;
    semilogy(forSave(:,10),forSave(:,14),'g');
    semilogy(forSave(:,11),forSave(:,15),'b');
    semilogy(forSave(:,12),forSave(:,16),'c');
    legend('npgs','npgs20','fistaBB','fistaL');
    keyboard

    system(['mv continuation.data traceLinGauss.data selectedTime.data timeVsA.data rmseVsA.data stepSizeLin.data varyMeasurement.data varyMeasurementTime.data skyline.data varyMeasurementTable.tex ' paperDir]);
    disp('done');
end

if(any(runList==903))
    filename = [mfilename '_003.mat'];
    load(filename);
    fprintf('PET Poisson example\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    a  = [0 0 0 0 0 -1];
    as = [2 2 2 1 1  1];

    K = 12;

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    npgscTime = mean(Cell.getField( npgsc(:,:,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,:,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    npgscCost = mean(Cell.getField( npgsc(:,:,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,:,1:K),'cost'),3);

    fbpRMSE   = mean(Cell.getField(   fbp(:,:,1:K),'RMSE'),3);
    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,:,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,:,1:K),'RMSE'),3);

    keyboard

    figure;
    loglog(count,   npgRMSE,'r-*'); hold on;
    loglog(count,   fbpRMSE,'b-o');
    loglog(count,spiralRMSE,'k-^');
    loglog(count,  npgcRMSE,'k*-.');
    loglog(count,  npgsRMSE,'c>-');
    loglog(count, npgscRMSE,'gs-');
    legend('npg','fbp','spiral','npgc','npgs','npgsc');

    figure;
    loglog(count,   npgTime,'r-*'); hold on;
    loglog(count,spiralTime,'k-^');
    loglog(count,  npgcTime,'k*-.');
    loglog(count,  npgsTime,'c>-');
    loglog(count, npgscTime,'gs-');
    legend('npg','spiral','npgc','npgs','npgsc');

    forSave=[npgTime, npgcTime, npgsTime, npgscTime, spiralTime,...
             npgCost, npgcCost, npgsCost, npgscCost, spiralCost,...
             npgRMSE, npgcRMSE, npgsRMSE, npgscRMSE, spiralRMSE,...
        fbpRMSE, count(:)];
    save('varyCntPET.data','forSave','-ascii');

    forSave=[]; t=0; mIdx=6;
    out=  npgc{mIdx,1,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=   npg{mIdx,1,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=spiral{mIdx,1,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=npgs{mIdx,1,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=npgsc{mIdx,1,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    mincost=reshape(forSave(:,[1,4,7]),[],1); 
    mincost=min(mincost(mincost~=0));
    idx=(forSave(:,1)~=0); forSave(idx,1)=(forSave(idx,1)-mincost);
    idx=(forSave(:,4)~=0); forSave(idx,4)=(forSave(idx,4)-mincost);
    idx=(forSave(:,7)~=0); forSave(idx,7)=(forSave(idx,7)-mincost);
    save('cost_itrPET.data','forSave','-ascii');

    figure; semilogy(forSave(:,3),forSave(:,1),'r'); hold on;
    semilogy(forSave(:,6),forSave(:,4),'g'); semilogy(forSave(:,9),forSave(:,7),'b');
    legend('npgc','npg','spiral');
    figure; semilogy(forSave(:,3),forSave(:,2),'r'); hold on;
    semilogy(forSave(:,6),forSave(:,5),'g'); semilogy(forSave(:,9),forSave(:,8),'b');
    legend('npgc','npg','spiral');

    nn=128;
    xtrue = read_zubal_emis('nx', nn, 'ny', nn);
    % attenuation map
    mumap = read_zubal_attn('nx', nn, 'ny', nn);
    imwrite(xtrue/max(xtrue(:)),'pet.png');
    imwrite(mumap/max(mumap(:)),'mumap.png');

    idx=5+6;
    fprintf('   NPG: %g%%\n',   npg{idx}.RMSE(end)*100);
    fprintf('  NPGc: %g%%\n',  npgc{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueAlpha)*100);
    fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
    fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
    fprintf('NPGscS: (%g%%, %g%%)\n',    npgsc_s.RMSE(end)*100,rmseTruncate(npgsc_s   )*100);
    img=npg{idx}.alpha; mask=npg{idx}.opt.mask;
    img=showImgMask(   npg{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPG_pet.png')
    img=showImgMask(  npgc{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGc_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGc_pet.png')
    img=showImgMask(spiral{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet.png')
    img=showImgMask(   fbp{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet.png')
    img=showImgMask(  npgs{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGs_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGs_pet.png')
    img=showImgMask( npgsc{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGsc_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGsc_pet.png')
    img=showImgMask( npgsc_s.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGscS_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGscS_pet.png')

    idx=4;
    fprintf('   NPG: %g%%\n',   npg{idx}.RMSE(end)*100);
    fprintf('  NPGc: %g%%\n',  npgc{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueAlpha)*100);
    fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
    fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
    img=npg{idx}.alpha; mask=npg{idx}.opt.mask;
    img=showImgMask(   npg{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPG_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPG_pet2.png')
    img=showImgMask(  npgc{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGc_pet2.png')
    img=showImgMask(spiral{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet2.png')
    img=showImgMask(   fbp{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet2.png')
    img=showImgMask(  npgs{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGs_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGs_pet2.png')
    img=showImgMask( npgsc{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGsc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGsc_pet2.png')
    
    keyboard
    system(['mv varyCntPET.data cost_itrPET.data ' paperDir]);
    close all;
end

if(any(runList==904))
    filename = [mfilename '_004.mat'];
    load(filename);

    keyboard

    prjFull = [60, 80, 100, 120, 180, 360];
    a=        [-6, -6,  -6,  -6,  -6,  -6];
    idx=2:2:7;
    K = 2;

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    npgscTime = mean(Cell.getField( npgsc(:,:,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,:,1:K),'time'),3);
    fpcasTime = mean(Cell.getField( fpcas(:,:,1:K),'time' ),3);
    fistaTime = mean(Cell.getField( fista(:,:,1:K),'time'),3);
    sparsaTime= mean(Cell.getField(sparsa(:,:,1:K),'time'),3);
    sparsnTime= mean(Cell.getField(sparsn(:,:,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    npgscCost = mean(Cell.getField( npgsc(:,:,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,:,1:K),'cost'),3);
    fpcasCost = mean(Cell.getField( fpcas(:,:,1:K),'f'   ),3);
    fistaCost = mean(Cell.getField( fista(:,:,1:K),'cost'),3);
    sparsaCost= mean(Cell.getField(sparsa(:,:,1:K),'cost'),3);
    sparsnCost= mean(Cell.getField(sparsn(:,:,1:K),'cost'),3);

    fbpRMSE   = mean(Cell.getField(   fbp(:,:,1:K),'RMSE'),3);
    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,:,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,:,1:K),'RMSE'),3);
    fpcasRMSE = mean(Cell.getField( fpcas(:,:,1:K),'RMSE'),3);
    fistaRMSE = mean(Cell.getField( fista(:,:,1:K),'RMSE'),3);
    sparsaRMSE= mean(Cell.getField(sparsa(:,:,1:K),'RMSE'),3);
    sparsnRMSE= mean(Cell.getField(sparsn(:,:,1:K),'RMSE'),3);

%   K=5;
%   npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
%   fpcasTime = mean(Cell.getField( fpcas(:,:,1:K),'time'),3);
%   npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
%   fpcasCost = mean(Cell.getField( fpcas(:,:,1:K),'f'),3);
%   npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
%   fpcasRMSE = mean(Cell.getField( fpcas(:,:,1:K),'RMSE'),3);
%   K=2;

    figure;
    semilogy(prjFull/2,   npgRMSE(:,1),'r-*'); hold on;
    semilogy(prjFull/2,  npgcRMSE(:,1),'c-p');
    semilogy(prjFull/2,  npgsRMSE(:,1),'k-s');
    semilogy(prjFull/2, npgscRMSE(:,1),'r-.');
    semilogy(prjFull/2, fpcasRMSE(:,1),'g-o');
    semilogy(prjFull/2,sparsaRMSE(:,1),'y-p');
    semilogy(prjFull/2,sparsnRMSE(:,1),'r-x');
    semilogy(prjFull/2,spiralRMSE(:,1),'k-.>');
    semilogy(prjFull/2,   fbpRMSE(:,1),'k-.s');
    semilogy(prjFull/2, fistaRMSE(:,1),'k-.o');
    legend('npg','npgc','npgs','npgsc','fpcas','sparsa','sparsn','spiral','fbp','fista');
    figure;
    semilogy(prjFull/2,   npgTime(:,1),'r-*'); hold on;
    semilogy(prjFull/2,  npgcTime(:,1),'c-p');
    semilogy(prjFull/2,  npgsTime(:,1),'k-s');
    semilogy(prjFull/2, npgscTime(:,1),'r-.');
    semilogy(prjFull/2, fpcasTime(:,1),'g-o');
    semilogy(prjFull/2,sparsaTime(:,1),'y-p');
    semilogy(prjFull/2,sparsnTime(:,1),'r-x');
    semilogy(prjFull/2,spiralTime(:,1),'k-.>');
    semilogy(prjFull/2, fistaTime(:,1),'k-.o');
    legend('npg','npgc','npgs','npgsc','fpcas','sparsa','sparsn','spiral','fista');

    snrIdx=[2;3;1;4];
    for i=3
        figure;
        semilogy(log10(snr(snrIdx)),   npgRMSE(i,snrIdx),'r-*'); hold on;
        semilogy(log10(snr(snrIdx)),  npgcRMSE(i,snrIdx),'c-p');
        semilogy(log10(snr(snrIdx)),  npgsRMSE(i,snrIdx),'k-s');
        semilogy(log10(snr(snrIdx)), fpcasRMSE(i,snrIdx),'g-o');
        semilogy(log10(snr(snrIdx)),sparsaRMSE(i,snrIdx),'y-p');
        semilogy(log10(snr(snrIdx)),sparsnRMSE(i,snrIdx),'r-x');
        semilogy(log10(snr(snrIdx)),spiralRMSE(i,snrIdx),'b-.>');
        semilogy(log10(snr(snrIdx)),   fbpRMSE(i,snrIdx),'k-.s');
        semilogy(log10(snr(snrIdx)), fistaRMSE(i,snrIdx),'k-.o');
        legend('npg','npgc','npgs','fpcas','sparsa','sparsn','spiral','fbp','fista');
        figure;
        semilogy(log10(snr(snrIdx)),   npgTime(i,snrIdx),'r-*'); hold on;
        semilogy(log10(snr(snrIdx)),  npgcTime(i,snrIdx),'c-p');
        semilogy(log10(snr(snrIdx)),  npgsTime(i,snrIdx),'k-s');
        semilogy(log10(snr(snrIdx)), fpcasTime(i,snrIdx),'g-o');
        semilogy(log10(snr(snrIdx)),sparsaTime(i,snrIdx),'y-p');
        semilogy(log10(snr(snrIdx)),sparsnTime(i,snrIdx),'r-x');
        semilogy(log10(snr(snrIdx)),spiralTime(i,snrIdx),'b-.>');
        semilogy(log10(snr(snrIdx)), fistaTime(i,snrIdx),'b-.o');
        legend('npg','npgc','npgs','fpcas','sparsa','sparsn','spiral','fista');
        keyboard;
    end

    forSave=prjFull(:);
    forSave=[forSave,   npgTime(:,1),   npgCost(:,1),   npgRMSE(:,1)];
    forSave=[forSave,  npgcTime(:,1),  npgcCost(:,1),  npgcRMSE(:,1)];
    forSave=[forSave,  npgsTime(:,1),  npgsCost(:,1),  npgsRMSE(:,1)];
    forSave=[forSave, fpcasTime(:,1), fpcasCost(:,1), fpcasRMSE(:,1)];
    forSave=[forSave,sparsaTime(:,1),sparsaCost(:,1),sparsaRMSE(:,1)];
    forSave=[forSave,sparsnTime(:,1),sparsnCost(:,1),sparsnRMSE(:,1)];
    forSave=[forSave,spiralTime(:,1),spiralCost(:,1),spiralRMSE(:,1)];
    forSave=[forSave, fistaTime(:,1), fistaCost(:,1), fistaRMSE(:,1)];
    forSave=[forSave, npgscTime(:,1), npgscCost(:,1), npgscRMSE(:,1)];
    forSave=[forSave,   fbpRMSE(:,1)];
    save('varyMeasurementPhantom.data','forSave','-ascii');

    forSave=snr(snrIdx); idx=3;
    forSave=[forSave;   npgTime(idx,snrIdx);   npgCost(idx,snrIdx);   npgRMSE(idx,snrIdx)];
    forSave=[forSave;  npgcTime(idx,snrIdx);  npgcCost(idx,snrIdx);  npgcRMSE(idx,snrIdx)];
    forSave=[forSave;  npgsTime(idx,snrIdx);  npgsCost(idx,snrIdx);  npgsRMSE(idx,snrIdx)];
    forSave=[forSave; fpcasTime(idx,snrIdx); fpcasCost(idx,snrIdx); fpcasRMSE(idx,snrIdx)];
    forSave=[forSave;sparsaTime(idx,snrIdx);sparsaCost(idx,snrIdx);sparsaRMSE(idx,snrIdx)];
    forSave=[forSave;sparsnTime(idx,snrIdx);sparsnCost(idx,snrIdx);sparsnRMSE(idx,snrIdx)];
    forSave=[forSave;spiralTime(idx,snrIdx);spiralCost(idx,snrIdx);spiralRMSE(idx,snrIdx)];
    forSave=[forSave; fistaTime(idx,snrIdx); fistaCost(idx,snrIdx); fistaRMSE(idx,snrIdx)];
    forSave=[forSave; npgscTime(idx,snrIdx); npgscCost(idx,snrIdx); npgscRMSE(idx,snrIdx)];
    forSave=[forSave;   fbpRMSE(idx,snrIdx)];
    forSave=forSave';
    save('varyNoisePhantom.data','forSave','-ascii');

    mIdx=3; as=1; forSave=[]; t=0;
    q=(1:max(find(sparsn12{mIdx,as}.cost(:)>=sparsn{mIdx,as}.cost(end))))';
    t=t+1; temp=      npg{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp; % 5th cal
    t=t+1; temp=     npgc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp; % 9th col
    t=t+1; temp= sparsn12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    q=(1:max(find(spiral12{mIdx,as}.cost(:)>=spiral{mIdx,as}.cost(end))))';
    t=t+1; temp= spiral12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp; % 17th col
    t=t+1; temp= spiral12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    q=(1:max(find(sparsa12{mIdx,as}.cost(:)>=sparsa{mIdx,as}.cost(end))))';
    t=t+1; temp= sparsa12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp; % 29th col
    t=t+1; temp= sparsa12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp=    fista{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp; % 33rd col
    t=t+1; temp=    fista{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    fista{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    fista{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp; % 37th col
    t=t+1; temp= sparsa12{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsa12{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp; % 41st col
    t=t+1; temp= spiral12{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= spiral12{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    save('traceLinGaussPhantom.data','forSave','-ascii');
    system(['mv varyMeasurementPhantom.data varyNoisePhantom.data traceLinGaussPhantom.data ' paperDir]);

    idx=3;
    fprintf('NPG: %g%%\n', npgc{idx}.RMSE(end)*100);
    fprintf('SpaRSA: %g%%\n', sparsn{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n', spiral{idx}.RMSE(end)*100);
    fprintf('NPGs: (%g%%, %g%%)\n',npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
    fprintf('SpaRSAs: (%g%%, %g%%)\n',sparsa{idx}.RMSE(end)*100,rmseTruncate(sparsa{idx},npgs{idx}.opt.trueAlpha)*100);
    fprintf('FBP: (%g%%, %g%%)\n',fbp{idx}.RMSE(end)*100,rmseTruncate(fbp{idx},npgs{idx}.opt.trueAlpha)*100);
    img=npgs{idx}.alpha; mask=npgs{idx}.opt.mask;
    img=showImgMask(npgc{idx}.alpha,mask); maxImg=max(img(:));
    figure; showImg(img,0); saveas(gcf,'NPG_phantom.eps','psc2');
    imwrite(img,'NPG_phantom.png')
    img=showImgMask(npgsc{idx}.alpha,mask); maxImg=max(img(:)); mask=npgsc{idx}.opt.mask;
    figure; showImg(img,0); saveas(gcf,'NPGs_phantom.eps','psc2');
    imwrite(img,'NPGs_phantom.png')
    img=showImgMask(fbp{idx}.alpha,mask); maxImg=max(img(:)); mask=npgs{idx}.opt.mask;
    figure; showImg(img,0); saveas(gcf,'FBP_phantom.eps','psc2');
    imwrite(img,'FBP_phantom.png')
    img=showImgMask(sparsn{idx}.alpha,mask); maxImg=max(img(:)); mask=npgs{idx}.opt.mask;
    figure; showImg(img,0); saveas(gcf,'SpaRSA_phantom.eps','psc2');
    imwrite(img,'SpaRSA_phantom.png')
    img=showImgMask(spiral{idx}.alpha,mask); maxImg=max(img(:)); mask=npgs{idx}.opt.mask;
    figure; showImg(img,0); saveas(gcf,'SPIRAL_phantom.eps','psc2');
    imwrite(img,'SPIRAL_phantom.png')
    img=showImgMask(sparsa{idx}.alpha,mask); maxImg=max(img(:)); mask=npgs{idx}.opt.mask;
    figure; showImg(img,0); saveas(gcf,'SpaRSAs_phantom.eps','psc2');
    imwrite(img,'SpaRSAs_phantom.png')

    pause;
    system(['mv *_phantom.png *_phantom.eps ' paperDir]);
    close all;

    return;

    for i=1:6
        a=npg{i};
        b=sparsn{i};
        c=npgc{i};
        d=spiral{i};
        mc=min([a.cost(:)' b.cost(:)' c.cost(:)' d.cost(:)']);
        figure(91); semilogy(a.time,a.cost-mc,'r'); hold on;
        semilogy(b.time,b.cost-mc,'b'); semilogy(c.time,c.cost-mc,'g');
        semilogy(d.time,d.cost-mc,'c');
        legend('NPG','SpaRSAp','NPGc','spiral');
        xlim([0,max(a.time(end),b.time(end))]);
        hold off;
        figure(92); semilogy(a.time,a.RMSE,'r'); hold on;
        semilogy(b.time,b.RMSE,'b');
        semilogy(c.time,c.RMSE,'g');
        semilogy(d.time,d.RMSE,'c');
        legend('NPG','SpaRSAp','NPGc','spiral');
        hold off;
        xlim([0,max(a.time(end),b.time(end))]);
        pause;
    end

    for i=1:6
        a=npgs{i};
        b=sparsa{i};
        c=npgsc{i};
        d=fpcas{i};
        mc=min([a.cost(:)' b.cost(:)' c.cost(:)' d.cost(:)']);
        figure(91); semilogy(a.time,a.cost-mc,'r'); hold on;
        semilogy(b.time,b.cost-mc,'b');
        semilogy(c.time,c.cost-mc,'g');
        semilogy(d.time,d.cost-mc,'k');
        legend('NPGs','SpaRSA','NPGsc','FPCas');
        hold off;
        figure(92); semilogy(a.time,a.RMSE,'r'); hold on;
        semilogy(b.time,b.RMSE,'b');
        semilogy(c.time,c.RMSE,'g');
        semilogy(d.time,d.RMSE,'k');
        legend('NPGs','SpaRSA','NPGsc','FPCas');
        hold off;
        pause;
    end

end

if(any(runList==905))
    filename = [mfilename '_005.mat']; load(filename);
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7]';
    u  =[1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3,1e-4,1e-4,1e-5];
    K=5;

    npgTime   = mean(Cell.getField(   npg(:,1,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,1,1:K),'time'),3);
    npgscTime = mean(Cell.getField( npgsc(:,1,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,1,1:K),'time'),3);
    fpcasTime = mean(Cell.getField( fpcas(:,1,1:K),'cpu' ),3);
    fistaTime = mean(Cell.getField( fista(:,1,1:K),'time'),3);
    sparsaTime= mean(Cell.getField(sparsa(:,1,1:K),'time'),3);
    sparsnTime= mean(Cell.getField(sparsn(:,1,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,1,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,1,1:K),'cost'),3);
    npgscCost = mean(Cell.getField( npgsc(:,1,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,1,1:K),'cost'),3);
    fpcasCost = mean(Cell.getField( fpcas(:,1,1:K),'f'   ),3);
    fistaCost = mean(Cell.getField( fista(:,1,1:K),'cost'),3);
    sparsaCost= mean(Cell.getField(sparsa(:,1,1:K),'cost'),3);
    sparsnCost= mean(Cell.getField(sparsn(:,1,1:K),'cost'),3);

    npgRMSE   = mean(Cell.getField(   npg(:,1,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,1,1:K),'RMSE'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,1,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,1,1:K),'reconerror'),3);
    fpcasRMSE = mean(Cell.getField( fpcas(:,1,1:K),'RMSE'),3);
    fistaRMSE = mean(Cell.getField( fista(:,1,1:K),'RMSE'),3);
    sparsaRMSE= mean(Cell.getField(sparsa(:,1,1:K),'RMSE'),3);
    sparsnRMSE= mean(Cell.getField(sparsn(:,1,1:K),'RMSE'),3);

    figure;
    loglog(snr,   npgRMSE,'r-*'); hold on;
    loglog(snr,  npgcRMSE,'c-p');
    loglog(snr, npgscRMSE,'k-s');
    loglog(snr,spiralRMSE,'k-^');
    loglog(snr, fpcasRMSE,'g-o');
    loglog(snr, fistaRMSE,'b-.');
    loglog(snr,sparsaRMSE,'y-p');
    loglog(snr,sparsnRMSE,'r-x');
    legend('npg','npgc','npgsc','spiral','fpcas','fista','sparsa','sparsn');
    figure;
    loglog(snr,   npgTime,'r-*'); hold on;
    loglog(snr,  npgcTime,'c-p');
    loglog(snr, npgscTime,'k-s');
    loglog(snr,spiralTime,'k-^');
    loglog(snr, fpcasTime,'g-o');
    loglog(snr, fistaTime,'b-.');
    loglog(snr,sparsaTime,'y-p');
    loglog(snr,sparsnTime,'r-x');
    legend('npg','npgc','npgsc','spiral','fpcas','fista','sparsa','sparsn');

    M=length(snr);
    str=        'SNR            ';              for i=1:M; str=sprintf('%s&%10d',str,snr(i)); end;
    str=sprintf('%s\\\\\\hline',str);
    str=sprintf('%s\\\\\nNPG            ', str);for i=1:M;str=sprintf('%s&%-10.4g',str,   npg{i}.cost(end));end;
    str=sprintf('%s\\\\\nNPG$_\\text{C}$ ',str);for i=1:M;str=sprintf('%s&%-10.4g',str,  npgc{i}.cost(end));end;
    str=sprintf('%s\\\\\nNPG$_\\text{S}$ ',str);for i=1:M;str=sprintf('%s&%-10.4g',str, npgsc{i}.cost(end));end;
    str=sprintf('%s\\\\\nSPIRAL         ', str);for i=1:M;str=sprintf('%s&%-10.4g',str,spiral{i}.cost(end));end;
    str=sprintf('%s\\\\\nFPC$_\\text{AS}$',str);for i=1:M;str=sprintf('%s&%-10.4g',str, fpcas{i}.cost(end));end;
    str=sprintf('%s\nFISTA          ', str);    for i=1:M;str=sprintf('%s&%-10.4g',str, fista{i}.cost(end));end;
    str=sprintf('%s\\\\\nSpaRSA         ', str);for i=1:M;str=sprintf('%s&%-10.4g',str,spiral{i}.cost(end));end;
    file=fopen('varySNRTable.tex','w'); fprintf(file,'%s',str); fclose(file);

      npgItr=[];   
     npgcItr=[];
    npgscItr=[];
   spiralItr=[];
    fpcasItr=[];
    fistaItr=[];
   sparsaItr=[];

    for i=1:K
        temp=   npg(:,1,i);    npgItr=[   npgItr,showResult(temp,2,'p'   )];
        temp=  npgc(:,1,i);   npgcItr=[  npgcItr,showResult(temp,2,'p'   )];
        temp= npgsc(:,1,i);  npgscItr=[ npgscItr,showResult(temp,2,'p'   )];
        temp=spiral(:,1,i); spiralItr=[spiralItr,showResult(temp,2,'p'   )];
        temp= fpcas(:,1,i);  fpcasItr=[ fpcasItr,showResult(temp,2,'itr' )];
        temp= fista(:,1,i);  fistaItr=[ fistaItr,showResult(temp,2,'p'   )];
        temp=sparsa(:,1,i); sparsaItr=[sparsaItr,showResult(temp,3,'RMSE')];
    end
    keyboard

    forSave=[    npgTime, npgcTime, npgscTime, spiralTime, fpcasTime, fistaTime, sparsaTime, sparsnTime, ...
                 npgCost, npgcCost, npgscCost, spiralCost, fpcasCost, fistaCost, sparsaCost, sparsnCost, ...
                 npgRMSE, npgcRMSE, npgscRMSE, spiralRMSE, fpcasRMSE, fistaRMSE, sparsaRMSE, sparsnRMSE, ...
                  snr(:)];
    save('varySNR.data','forSave','-ascii');

    system(['mv varySNRTable.tex varySNR.data ' paperDir]);
end

if(any(runList==906))
    filename = [mfilename '_006.mat']; load(filename);
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    aa =10.^(3:-1:-10);

    K=5;
    for k=1:K
        for i=1:length(m)
             npgContRMSE(:,i,k)  = [ npgFull{i,k}.contRMSE(:); npgFull{i,k}.RMSE(end)];
            npgsContRMSE(:,i,k)  = [npgsFull{i,k}.contRMSE(:);npgsFull{i,k}.RMSE(end)];
        end
    end
    npgContRMSE=mean(npgContRMSE,3);
    npgsContRMSE=mean(npgsContRMSE,3);

    figure;
    for i=1:length(m)
        loglog(aa,npgContRMSE(:,i),'r'); hold on; loglog(aa,npgsContRMSE(:,i),'b');
        aaa(i)=min(npgContRMSE(:,i));
        bbb(i)=min(npgsContRMSE(:,i));
    end
    figure; semilogy(aaa); hold on; semilogy(bbb);

    keyboard

    K=5;
    fprintf('Poisson example\n');

    npgTime    = mean(Cell.getField(    npg(:,:,1:K),'time'),3);
    npgcTime   = mean(Cell.getField(   npgc(:,:,1:K),'time'),3);
    npgsTime   = mean(Cell.getField(   npgs(:,:,1:K),'time'),3);
    npgscTime  = mean(Cell.getField(  npgsc(:,:,1:K),'time'),3);
    spiralTime = mean(Cell.getField( spiral(:,:,1:K),'time'),3);
    spiral8Time= mean(Cell.getField(spiral8(:,:,1:K),'time'),3);

    npgCost    = mean(Cell.getField(    npg(:,:,1:K),'cost'),3);
    npgcCost   = mean(Cell.getField(   npgc(:,:,1:K),'cost'),3);
    npgsCost   = mean(Cell.getField(   npgs(:,:,1:K),'cost'),3);
    npgscCost  = mean(Cell.getField(  npgsc(:,:,1:K),'cost'),3);
    spiralCost = mean(Cell.getField( spiral(:,:,1:K),'cost'),3);
    spiral8Cost= mean(Cell.getField(spiral8(:,:,1:K),'cost'),3);

    npgRMSE    = mean(Cell.getField(    npg(:,:,1:K),'RMSE'),3);
    npgcRMSE   = mean(Cell.getField(   npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE   = mean(Cell.getField(   npgs(:,:,1:K),'RMSE'),3);
    npgscRMSE  = mean(Cell.getField(  npgsc(:,:,1:K),'RMSE'),3);
    spiralRMSE = mean(Cell.getField( spiral(:,:,1:K),'RMSE'),3);
    spiral8RMSE= mean(Cell.getField(spiral8(:,:,1:K),'RMSE'),3);

    aIdx=4;
    figure;
    semilogy(m,    npgRMSE(:,aIdx),'r-*'); hold on;
    semilogy(m,   npgsRMSE(:,aIdx),'c-p');
    semilogy(m, spiralRMSE(:,aIdx),'k-^');
    semilogy(m,   npgcRMSE(:,aIdx),'k*-.');
    semilogy(m,  npgscRMSE(:,aIdx),'bs-.');
    semilogy(m,spiral8RMSE(:,aIdx),'go-.');
    legend('npg','npgs','spiral','npgc','npgsc','spiral8');

    figure;
    semilogy(m,    npgTime(:,aIdx),'r-*'); hold on;
    semilogy(m,   npgsTime(:,aIdx),'c-p');
    semilogy(m, spiralTime(:,aIdx),'k-^');
    semilogy(m,   npgcTime(:,aIdx),'k*-.');
    semilogy(m,  npgscTime(:,aIdx),'bs-.');
    semilogy(m,spiral8Time(:,aIdx),'go-.');
    legend('npg','npgs','spiral','npgc','npgsc','spiral8');

    forSave=[npgTime(:,aIdx), npgsTime(:,aIdx), npgcTime(:,aIdx), npgscTime(:,aIdx), spiralTime(:,aIdx), spiral8Time(:,aIdx),...
             npgCost(:,aIdx), npgsCost(:,aIdx), npgcCost(:,aIdx), npgscCost(:,aIdx), spiralCost(:,aIdx), spiral8Cost(:,aIdx),...
             npgRMSE(:,aIdx), npgsRMSE(:,aIdx), npgcRMSE(:,aIdx), npgscRMSE(:,aIdx), spiralRMSE(:,aIdx), spiral8RMSE(:,aIdx),...
        m(:)];
    save('varyMeasurementPoisson.data','forSave','-ascii');

    temp = 2;
    signal=npg{1}.opt.trueAlpha;
    signal=[signal,    npg{temp,3,1}.alpha];
    signal=[signal,   npgc{temp,3,1}.alpha];
    signal=[signal, spiral{temp,4,1}.alpha];
    signal=[signal,   npgs{temp,3,1}.alpha];
    signal=[signal,  npgsc{temp,3,1}.alpha];
    save('skylinePoisson.data','signal','-ascii');
    fprintf('   NPG: %g%%\n',           npg{temp,3,1}.RMSE(end)*100);
    fprintf('  NPGc: %g%%\n',          npgc{temp,3,1}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',        spiral{temp,4,1}.RMSE(end)*100);
    fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{temp,3,1}.RMSE(end)*100,rmseTruncate( npgs{temp,3,1})*100);
    fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{temp,3,1}.RMSE(end)*100,rmseTruncate(npgsc{temp,3,1})*100);
    system(['mv skylinePoisson.data ' paperDir]);

    forSave=[]; t=0; mIdx=2;
    out=  npgc{mIdx,aIdx,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=   npg{mIdx,aIdx,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=spiral8{mIdx,aIdx,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=   npgs{mIdx,aIdx,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=  npgsc{mIdx,aIdx,1};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=spiral8{mIdx,aIdx,1};
    q=(1:max(find(spiral8{mIdx,aIdx}.cost(:)>=spiral{mIdx,aIdx}.cost(end))))';
    t=t+1; forSave(1:length(out.cost(q)),t)=out.cost(q);
    t=t+1; forSave(1:length(out.RMSE(q)),t)=out.RMSE(q);
    t=t+1; forSave(1:length(out.time(q)),t)=out.time(q);

    save('cost_itr.data','forSave','-ascii');
    mincost=reshape(forSave(:,[1,4,7,10,13,16]),[],1); 
    mincost=min(mincost(mincost~=0));
    idx=(forSave(:, 1)~=0); forSave(idx, 1)=(forSave(idx, 1)-mincost);
    idx=(forSave(:, 4)~=0); forSave(idx, 4)=(forSave(idx, 4)-mincost);
    idx=(forSave(:, 7)~=0); forSave(idx, 7)=(forSave(idx, 7)-mincost);
    idx=(forSave(:,10)~=0); forSave(idx,10)=(forSave(idx,10)-mincost);
    idx=(forSave(:,13)~=0); forSave(idx,13)=(forSave(idx,13)-mincost);

    figure; semilogy(forSave(:,1),'r'); hold on;
    semilogy(forSave(:,4),'g'); semilogy(forSave(:,7),'b');
    semilogy(forSave(:,10),'k'); semilogy(forSave(:,13),'c');
    legend('npgc','npg','spiral8','npgs','npgsc');
    figure; loglog(forSave(:,3),forSave(:,1),'r'); hold on;
    semilogy(forSave(:,6),forSave(:,4),'g'); semilogy(forSave(:,9),forSave(:,7),'b');
    semilogy(forSave(:,12),forSave(:,10),'k'); semilogy(forSave(:,15),forSave(:,13),'c');
    legend('npgc','npg','spiral8','npgs','npgsc');
    figure; loglog(forSave(:,3),forSave(:,2),'r'); hold on;
    semilogy(forSave(:,6),forSave(:,5),'g'); semilogy(forSave(:,9),forSave(:,8),'b');
    semilogy(forSave(:,12),forSave(:,11),'k'); semilogy(forSave(:,15),forSave(:,14),'c');
    legend('npgc','npg','spiral8','npgs','npgsc');

    keyboard

    printf('\nstart to vary a to see the performance\n');

    % 3 groups are enough
    idx=[1 2 6];
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
    as=1:5;
    forSave=[]; forTime=[];
    for mIdx=idx
        figure(900);
        semilogy(log10(a(mIdx))+as-3,    npgRMSE(mIdx,as),'r-*'); hold on;
        semilogy(log10(a(mIdx))+as-3,   npgcRMSE(mIdx,as),'r.-');
        semilogy(log10(a(mIdx))+as-3, spiralRMSE(mIdx,as),'r-^');
        semilogy(log10(a(mIdx))+as-3,spiral8RMSE(mIdx,as),'g-o');
        semilogy(log10(a(mIdx))+as-3,   npgsRMSE(mIdx,as),'k-s');
        semilogy(log10(a(mIdx))+as-3,  npgscRMSE(mIdx,as),'b-p');
        legend('npg','npgc','spiral','spiral8','npgs','npgsc');

        forSave=[forSave log10(a(mIdx))+as(:)-3];
        forSave=[forSave reshape(    npgRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(   npgcRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape( spiralRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(spiral8RMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(   npgsRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(  npgscRMSE(mIdx,as),[],1)];

        figure(901);
        semilogy(log10(a(mIdx))+as-3,    npgTime(mIdx,as),'r-*'); hold on;
        semilogy(log10(a(mIdx))+as-3,   npgcTime(mIdx,as),'r.-');
        semilogy(log10(a(mIdx))+as-3, spiralTime(mIdx,as),'r-^');
        semilogy(log10(a(mIdx))+as-3,spiral8Time(mIdx,as),'g-o');
        semilogy(log10(a(mIdx))+as-3,   npgsTime(mIdx,as),'k-s');
        semilogy(log10(a(mIdx))+as-3,  npgscTime(mIdx,as),'b-p');
        legend('npg','npgc','spiral','spiral8','npgs','npgsc');
        title(sprintf('mIdx=%d',mIdx));

        forTime=[forTime log10(a(mIdx))+as(:)-3];
        forTime=[forTime reshape(    npgTime(mIdx,as),[],1)];
        forTime=[forTime reshape(   npgcTime(mIdx,as),[],1)];
        forTime=[forTime reshape( spiralTime(mIdx,as),[],1)];
        forTime=[forTime reshape(spiral8Time(mIdx,as),[],1)];
        forTime=[forTime reshape(   npgsTime(mIdx,as),[],1)];
        forTime=[forTime reshape(  npgscTime(mIdx,as),[],1)];
    end
    save('rmseVsAPoisson.data','forSave','-ascii');
    save('timeVsAPoisson.data','forTime','-ascii');

    system(['mv varyMeasurementPoisson.data cost_itr.data ' paperDir]);
end

if(any(runList==907))
    filename = [mfilename '_007.mat']; load(filename);

    m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200
    for k=1:4
        for i=1:length(m)
            npgRMSE(i,k) = min(npgFull{i,k}.contRMSE);
            npgsRMSE(i,k) = min(npgsFull{i,k}.contRMSE);

            gnetRMSE(i,k)=mean(gnet{i,k}.RMSE(:));

            npg0RMSE(i,k) = min(npgFull_knownI0{i,k}.contRMSE);
            npgs0RMSE(i,k) = min(npgsFull_knownI0{i,k}.contRMSE);

            temp= npglwFull{i,k}; if(~isempty(temp))  npglwRMSE(i,k)=min(temp.contRMSE); else  npglwRMSE(i,k)=0; end
            temp=npgslwFull{i,k}; if(~isempty(temp)) npgslwRMSE(i,k)=min(temp.contRMSE); else npgslwRMSE(i,k)=0; end
            % temp=  npglFull{i,k}; if(~isempty(temp))   npglRMSE(i,k)=min(temp.contRMSE); else   npglRMSE(i,k)=0; end
            % temp= npgslFull{i,k}; if(~isempty(temp))  npgslRMSE(i,k)=min(temp.contRMSE); else  npgslRMSE(i,k)=0; end
        end

        figure;
        semilogy(m,   npgRMSE(:,k),'r-*'); hold on;
        plot(    m,  npgsRMSE(:,k),'b-*');
        plot(    m,  npg0RMSE(:,k),'r.-');
        plot(    m, npgs0RMSE(:,k),'b.-');
        plot(    m, npglwRMSE(:,k),'rs-');
        plot(    m,npgslwRMSE(:,k),'bs-');
        % plot(    m,  npglRMSE,'rp-');
        % plot(    m, npgslRMSE,'bp-');
        legend('npg','npgs','npg0','npgs0','npglw','npwslw');
    end
end

if(any(runList==909))
    filename = [mfilename '_009.mat']; load(filename);

    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    fprintf('Poisson Log link example with glass beads\n');

    npgTime   = Cell.getField(   npg,'time');
    npgsTime  = Cell.getField(  npgs,'time');
    npgCost   = Cell.getField(   npg,'cost');
    npgsCost  = Cell.getField(  npgs,'cost');
    npgRMSE   = Cell.getField(   npg,'RMSE');
    npgsRMSE  = Cell.getField(  npgs,'RMSE');
    fbpRMSE   = Cell.getField(   fbp,'RMSE');

    for i=1:length(npgsRMSE(:,1))
        for j=1:length(npgsRMSE(1,:))
            npgsTranRMSE(i,j) = rmseTruncate(npgs{i,j});
        end
    end

    [r,c1]=find(   npgRMSE==repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c3]=find(  npgsTranRMSE==repmat(min(  npgsTranRMSE,[],2),1,5)); [r,idx3]=sort(r);
    disp([c1(idx1) ,c2(idx2), c3(idx3)]);
    idx1=(c1(idx1)-1)*6+(1:6)';
    idx2=(c2(idx2)-1)*6+(1:6)';
    idx3=(c3(idx3)-1)*6+(1:6)';
    idx2=idx1;

    figure;
    semilogy(prjFull,   npgRMSE(idx1),'r-*'); hold on;
    plot(prjFull,  npgsRMSE(idx2),'c-p');
    plot(prjFull,  npgsTranRMSE(idx3),'b.-');
    plot(prjFull, fbpRMSE(:),'g-s');
    legend('npg','npgs','npgsTran','fbp');

    figure;
    plot(prjFull,   npgTime(idx1),'r-*'); hold on;
    loglog(prjFull,  npgsTime(idx2),'c-p');
    %loglog(prjFull,   fbpTime(:,idx3),'g-s');
    legend('npg','npgs');

    for i=1:length(wnpgsFull)
        wnpgRMSE(i)=min(wnpgFull{i}.contRMSE);
        wnpgsRMSE(i)=min(wnpgsFull{i}.contRMSE);
    end

    forSave=[];
    forSave=[forSave,     npgRMSE(idx1)];
    forSave=[forSave,    npgsRMSE(idx2)];
    forSave=[forSave,     fbpRMSE(:)];
    forSave=[forSave,npgsTranRMSE(idx3)];

    forSave=[forSave,     npgTime(idx1)];
    forSave=[forSave,    npgsTime(idx2)];

    % forSave=[forSave,    npgCost(:,idx1)];
    % forSave=[forSave,   npgsCost(:,idx2)];
    % forSave=[forSave,    fbpCost(:,idx3)];

    forSave=[forSave,  prjFull(:)];
    forSave=[forSave,    wnpgRMSE(:)];
    forSave=[forSave,   wnpgsRMSE(:)];
    save('varyPrjGlassBead.data','forSave','-ascii');

    idx=5;
    img=showImgMask(npg {idx1(idx)}.alpha,npg{idx1(idx)}.opt.mask);
    maxImg=max(img(:));
    maxImg=1.2;
    figure; showImg(img);
    imwrite(img/maxImg,'NPGgb.png','png');
    img=showImgMask(npgs{idx2(idx)}.alpha,npg{idx1(idx)}.opt.mask);
    imwrite(img/maxImg,'NPGSgb.png','png'); figure; showImg(img,0,maxImg);
    img=showImgMask(fbp {idx}.alpha,npg{idx1(idx)}.opt.mask);
    imwrite(img/maxImg,'FBPgb.png','png'); figure; showImg(img,0,maxImg);
    img=showImgMask(opt.trueAlpha,npg{idx1(idx)}.opt.mask);
    imwrite(img/maxImg,'glassbeads.png','png');

    disp([npg{idx1(idx)}.RMSE(end), npgs{idx2(idx)}.RMSE(end), fbp{idx}.RMSE]);
    trueAlpha=npg{idx1(idx)}.opt.trueAlpha;
    disp([rmseTruncate(npg{idx1(idx)},trueAlpha), rmseTruncate(npgs{idx2(idx)},trueAlpha), rmseTruncate(fbp{idx},trueAlpha)]);

    system(['mv varyPrjGlassBead.data NPGgb.png NPGSgb.png FBPgb.png glassbeads.png ' paperDir]);
    keyboard
end

end

function e=gEle(x,i); e=x(i); end
function rmse=rmseTruncate(x,trueAlpha)
    if(nargin<2) trueAlpha=x.opt.trueAlpha; end;
    alpha=x.alpha; alpha(alpha<0)=0;
    rmse=sqrNorm(alpha-trueAlpha)/sqrNorm(trueAlpha);
end

