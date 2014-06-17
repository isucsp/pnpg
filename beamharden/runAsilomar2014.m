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
    opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1; %opt.matrixType='nonneg';
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    a=[1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5];
    snr=[inf 1e6 1e5 2 5 10 100 1e3 inf];
    for k=1:5
        for i=1:length(m)
            opt.m=m(i); opt.snr=1e6;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u(i) = a(i)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=5:-1:1
                fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
                opt.u = u(i)*10^(j-2);

                opt.continuation=false;
                opt.alphaStep='FISTA_ADMM_NNL1';
                %if(j<4) initSig = npg001{i,j+1,k}.alpha; end
                npg001{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npg001','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgC001{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgC001','-append');

                opt.continuation=false;
                opt.alphaStep='FISTA_L1';
                %if(j<4) initSig = FISTA001{i,j+1,k}.alpha; end
                FISTA001{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'FISTA001','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_L1';
                %if(j<4) initSig = FISTA001{i,j+1,k}.alpha; end
                FISTAC001{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'FISTAC001','-append');

                %if(j<4) initSig = fpcas001{i,j+1,k}.alpha; end
                A = @(xx) conf.Phi(conf.Psi(xx));
                At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                % option.gtol = 1e-3; option.gtol_scale_x = 1e-6; option.sub_mxitr = 10;
                [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas001{i,j,k}=out; fpcas001{i,j,k}.alpha = conf.Psi(s);
                fpcas001{i,j,k}.fVal(1)=0.5*norm(conf.Phi(fpcas001{i,j,k}.alpha)-conf.y)^2;
                fpcas001{i,j,k}.fVal(2)=norm(fpcas001{i,j,k}.alpha.*(fpcas001{i,j,k}.alpha<0))^2;
                fpcas001{i,j,k}.fVal(3)=sum(abs(conf.Psit(fpcas001{i,j,k}.alpha)));
                fpcas001{i,j,k}.opt = opt; alphaHat=fpcas001{i,j,k}.alpha;
                fpcas001{i,j,k}.RMSE=(norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha)).^2;
                fprintf('fpcas RMSE=%g\n',fpcas001{i,j,k}.RMSE);
                save(filename,'fpcas001','-append');
            end
        end
    end
end

% vary the number of the measurements and get the corresponding good u for 
% each, Use the matrix with positive elements as used in poisson example
if(any(runList==011))
    load(filename,'*011');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5;
    opt.thresh=1e-6;
    opt.matrixType='nonneg';
    m=[200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u = 10.^[1 0 -1 -2 -3 -4 -5 -6 -7];
    snr=[inf 1e6 1e5 2 5 10 100 1e3 inf];
    for i=4:4
        opt.m=m(i); opt.snr=inf;
        opt=loadLinear(conf,opt);
        for j=2:5
            fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
            opt.u = u(j);
            initSig = conf.Phit(conf.y)*0;

            opt.continuation=true;
            opt.alphaStep='FISTA_ADMM_NNL1';
            npgC011{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgC011','-append');

            opt.continuation=true;
            opt.alphaStep='FISTA_L1';
            FISTAC011{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'FISTAC011','-append');

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
            option.mxitr=opt.maxItr;
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            fpcas011{i,j}=out; fpcas011{i,j}.alpha = conf.Psi(s);
            fpcas011{i,j}.fVal(1)=0.5*norm(conf.Phi(fpcas011{i,j}.alpha)-conf.y)^2;
            fpcas011{i,j}.fVal(2)=norm(fpcas011{i,j}.alpha.*(fpcas011{i,j}.alpha<0))^2;
            fpcas011{i,j}.fVal(3)=sum(abs(conf.Psit(fpcas011{i,j}.alpha)));
            fpcas011{i,j}.opt = opt; alphaHat=fpcas011{i,j}.alpha;
            fpcas011{i,j}.RMSE=norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha);
            fprintf('fpcas RMSE=%g\n',fpcas011{i,j}.RMSE);
            save(filename,'fpcas011','-append');
        end
    end
end

% vary the inter loop criteria
if(any(runList==021))
    load(filename,'*021');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
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
        subplot(4,2,6); plot(npg021{1,j}.time,npg021{1,j}.time,style{t}); hold on;
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
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
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
        subplot(4,2,6); plot(npg031{1,j}.time,npg031{1,j}.time,style{t}); hold on;
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
    load(filename,'*002');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT();
    opt.maxItr=1e5; opt.thresh=1e-6; opt.debugLevel=0;
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-4,1e-4,1e-4,1e-4,1e-4,1e-4,1e-5,1e-5,1e-5];
    for k=1:2
        for i=1:length(m)
            opt.m=m(i); opt.snr=inf;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = (1e-5)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=1:5
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
                opt.u = u*10^(j-2);

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgC002{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgC002','-append');

                opt.continuation=false;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npg002{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npg002','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_L1';
                FISTAC002{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'FISTAC002','-append');

                opt.continuation=false;
                opt.alphaStep='FISTA_L1';
                FISTA002{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'FISTA002','-append');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral002{i,j,k}=out;
                spiral002{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral002{i,j,k}.alpha)-conf.y);
                spiral002{i,j,k}.fVal(2)=sqrNorm(spiral002{i,j,k}.alpha.*(spiral002{i,j,k}.alpha<0));
                spiral002{i,j,k}.fVal(3)=pNorm(conf.Psit(spiral002{i,j,k}.alpha),1);
                save(filename,'spiral002','-append');

                A = @(xx) conf.Phi(conf.Psi(xx));
                At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas002{i,j,k}=out; fpcas002{i,j,k}.alpha = conf.Psi(s);
                fpcas002{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(fpcas002{i,j,k}.alpha)-conf.y);
                fpcas002{i,j,k}.fVal(2)=sqrNorm(fpcas002{i,j,k}.alpha.*(fpcas002{i,j,k}.alpha<0));
                fpcas002{i,j,k}.fVal(3)=pNorm(conf.Psit(fpcas002{i,j,k}.alpha),1);
                fpcas002{i,j,k}.opt = opt; alphaHat=fpcas002{i,j,k}.alpha;
                fpcas002{i,j,k}.RMSE=sqrNorm(alphaHat-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
                fprintf('fpcas RMSE=%g\n',fpcas002{i,j,k}.RMSE);
                save(filename,'fpcas002','-append');
            end
        end
    end
end

if(any(runList==902))
    load(filename,'*002'); j=1;
    opt.a=[];
    opt=loadLinear(ConfigCT(),opt);
    signal=opt.trueAlpha;
    save('skyline.data','signal','-ascii');

    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    K = 2;

    npgTime=      0;
    npgCTime=     0;
    fistaCTime=   0;
    fistaTime=    0;
    spiralTime=   0;
    fpcasTime=    0;
    npgCost=      0;
    npgCCost=     0;
    fistaCCost=   0;
    fistaCost=    0;
    spiralCost=   0;
    fpcasCost=    0;
    npgRMSE=      0;
    npgCRMSE=     0;
    fistaCRMSE=   0;
    fistaRMSE=    0;
    spiralRMSE=   0;
    fpcasRMSE=    0;

    for k=1:K
        npgTime=        npgTime+1/K*showResult(npg002(:,:,k),2,'time');
        npgCTime=        npgCTime+1/K*showResult(npgC002(:,:,k),2,'time');
        fistaCTime=        fistaCTime+1/K*showResult(FISTAC002(:,:,k),2,'time');
        fistaTime=        fistaTime+1/K*showResult(FISTA002(:,:,k),2,'time');
        spiralTime=        spiralTime+1/K*showResult(spiral002(:,:,k),2,'time');
        fpcasTime=        fpcasTime+1/K*showResult(fpcas002(:,:,k),2,'cpu');

        npgCost=        npgCost+1/K*showResult(npg002(:,:,k),2,'cost');
        npgCCost=        npgCCost+1/K*showResult(npgC002(:,:,k),2,'cost');
        fistaCCost=        fistaCCost+1/K*showResult(FISTAC002(:,:,k),2,'cost');
        fistaCost=        fistaCost+1/K*showResult(FISTA002(:,:,k),2,'cost');
        spiralCost=        spiralCost+1/K*showResult(spiral002(:,:,k),2,'cost');
        fpcasCost=        fpcasCost+1/K*showResult(fpcas002(:,:,k),2,'f');

        npgRMSE=        npgRMSE+1/K*showResult(npg002(:,:,k),2,'RMSE');
        npgCRMSE=        npgCRMSE+1/K*showResult(npgC002(:,:,k),2,'RMSE');
        fistaCRMSE=        fistaCRMSE+1/K*showResult(FISTAC002(:,:,k),2,'RMSE');
        fistaRMSE=        fistaRMSE+1/K*showResult(FISTA002(:,:,k),2,'RMSE');
        spiralRMSE=        spiralRMSE+1/K*showResult(spiral002(:,:,k),2,'reconerror');
        fpcasRMSE=        fpcasRMSE+1/K*showResult(fpcas002(:,:,k),2,'RMSE');
    end

    % npgTime(:,1:2)=[];
    % npgCTime(:,1:2)=[];
    % fistaCTime(:,1:2)=[];
    % spiralTime(:,1:2)=[];
    % fpcasTime(:,1:2)=[];
    % npgCost(:,1:2)=[];
    % npgCCost(:,1:2)=[];
    % fistaCCost(:,1:2)=[];
    % spiralCost(:,1:2)=[];
    % fpcasCost(:,1:2)=[];
    % npgRMSE(:,1:2)=[];
    % npgCRMSE(:,1:2)=[];
    % fistaCRMSE(:,1:2)=[];
    % spiralRMSE(:,1:2)=[];
    % fpcasRMSE(:,1:2)=[];
    % fpcasRMSE_(:,1:2)=[];
    % fistaCRMSE_(:,1:2)=[];

    [r,c1]=find(npgRMSE  == repmat(min(npgRMSE   ,[],2),1,5)); [r,idx1]=sort(r); [r,c1(idx1)]
    [r,c2]=find(npgCRMSE == repmat(min(npgCRMSE  ,[],2),1,5)); [r,idx2]=sort(r); [r,c2(idx2)]
    [r,c3]=find(fistaCRMSE==repmat(min(fistaCRMSE,[],2),1,5)); [r,idx3]=sort(r); [r,c3(idx3)]
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r); [r,c4(idx4)]
    [r,c5]=find(fpcasRMSE== repmat(min(fpcasRMSE ,[],2),1,5)); [r,idx5]=sort(r); [r,c5(idx5)]
    [r,c6]=find(fistaRMSE== repmat(min(fistaRMSE ,[],2),1,5)); [r,idx6]=sort(r); [r,c6(idx6)]

    figure;
    semilogy(m,npgRMSE   ((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,npgCRMSE  ((c1(idx1)-1)*9+(1:9)'),'c-p');
    semilogy(m,fistaCRMSE((c1(idx1)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c1(idx1)-1)*9+(1:9)'),'k-^');
    semilogy(m,fpcasRMSE ((c1(idx1)-1)*9+(1:9)'),'g-o');
    semilogy(m,fistaRMSE ((c1(idx1)-1)*9+(1:9)'),'b-.');
    figure;
    semilogy(m,npgTime   ((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,npgCTime  ((c1(idx1)-1)*9+(1:9)'),'c-p');
    semilogy(m,fistaCTime((c1(idx1)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c1(idx1)-1)*9+(1:9)'),'k-^');
    semilogy(m,fpcasTime ((c1(idx1)-1)*9+(1:9)'),'g-o');
    semilogy(m,fistaTime ((c1(idx1)-1)*9+(1:9)'),'b-.');
    
    figure;
    semilogy(m,npgRMSE   ((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,npgCRMSE  ((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,fistaCRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m,fpcasRMSE ((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m,fistaRMSE ((c6(idx6)-1)*9+(1:9)'),'b-.');
    figure;
    semilogy(m,npgTime   ((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,npgCTime  ((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,fistaCTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m,fpcasTime ((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m,fistaTime ((c6(idx6)-1)*9+(1:9)'),'b-.');
    keyboard
    figure;
    semilogy(m,npgRMSE   (:,2),'r-*'); hold on;
    semilogy(m,npgCRMSE  (:,2),'c-p');
    semilogy(m,fistaCRMSE(:,2),'k-s');
    semilogy(m,spiralRMSE(:,2),'k-^');
    semilogy(m,fpcasRMSE (:,2),'g-o');
    semilogy(m,fistaRMSE (:,2),'b-.');
    

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
    forSave=[forSave, sum(fistaCRMSE_,2)./sum(fistaCRMSE_>0,2)];
    forSave=[forSave, sum(fpcasRMSE_,2)./sum(fpcasRMSE_>0,2)];
    save('varyMeasurement.data','forSave','-ascii');
end

% vary the number of measurements, with continuation, nonneg Phi
if(any(runList==012))
    load(filename,'*012');
    s = RandStream.create('mt19937ar','seed',0);
    RandStream.setGlobalStream(s);
    conf=ConfigCT(); opt.debugLevel=0;
    opt.maxItr=1e5; opt.thresh=1e-6; opt.matrixType='gaussian';
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    for k=1:2
        for i=1:length(m)
            opt.m=m(i); opt.snr=1e6;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = (1e-5)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=1:5
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
                opt.u = u*10^(j-2);

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgC012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgC012','-append');

                opt.continuation=false;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npg012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npg012','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_L1';
                FISTAC012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'FISTAC012','-append');

                opt.continuation=false;
                opt.alphaStep='FISTA_L1';
                FISTA012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'FISTA012','-append');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral012{i,j,k}=out;
                spiral012{i,j,k}.fVal(1)=0.5*norm(conf.Phi(spiral012{i,j,k}.alpha)-conf.y)^2;
                spiral012{i,j,k}.fVal(2)=norm(spiral012{i,j,k}.alpha.*(spiral012{i,j,k}.alpha<0))^2;
                spiral012{i,j,k}.fVal(3)=sum(abs(conf.Psit(spiral012{i,j,k}.alpha)));
                save(filename,'spiral012','-append');

                A = @(xx) conf.Phi(conf.Psi(xx));
                At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas012{i,j,k}=out; fpcas012{i,j,k}.alpha = conf.Psi(s);
                fpcas012{i,j,k}.fVal(1)=0.5*norm(conf.Phi(fpcas012{i,j,k}.alpha)-conf.y)^2;
                fpcas012{i,j,k}.fVal(2)=norm(fpcas012{i,j,k}.alpha.*(fpcas012{i,j,k}.alpha<0))^2;
                fpcas012{i,j,k}.fVal(3)=sum(abs(conf.Psit(fpcas012{i,j,k}.alpha)));
                fpcas012{i,j,k}.opt = opt; alphaHat=fpcas012{i,j,k}.alpha;
                fpcas012{i,j,k}.RMSE=(norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha))^2;
                fprintf('fpcas RMSE=%g\n',fpcas012{i,j,k}.RMSE);
                save(filename,'fpcas012','-append');
            end
        end
    end
end

if(any(runList==912))
    load(filename,'*012'); j=1;
    opt.a=[];
    opt=loadLinear(ConfigCT(),opt);
    signal=opt.trueAlpha;
    save('skyline.data','signal','-ascii');

    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    K = 1;

    npgTime=      0;
    npgCTime=     0;
    fistaCTime=   0;
    fistaTime=    0;
    spiralTime=   0;
    fpcasTime=    0;
    npgCost=      0;
    npgCCost=     0;
    fistaCCost=   0;
    fistaCost=    0;
    spiralCost=   0;
    fpcasCost=    0;
    npgRMSE=      0;
    npgCRMSE=     0;
    fistaCRMSE=   0;
    fistaRMSE=    0;
    spiralRMSE=   0;
    fpcasRMSE=    0;

    for k=1:K
        npgTime=        npgTime+1/K*showResult(npg012(:,:,k),2,'time');
        npgCTime=        npgCTime+1/K*showResult(npgC012(:,:,k),2,'time');
        fistaCTime=        fistaCTime+1/K*showResult(FISTAC012(:,:,k),2,'time');
        fistaTime=        fistaTime+1/K*showResult(FISTA012(:,:,k),2,'time');
        spiralTime=        spiralTime+1/K*showResult(spiral012(:,:,k),2,'time');
        fpcasTime=        fpcasTime+1/K*showResult(fpcas012(:,:,k),2,'cpu');

        npgCost=        npgCost+1/K*showResult(npg012(:,:,k),2,'cost');
        npgCCost=        npgCCost+1/K*showResult(npgC012(:,:,k),2,'cost');
        fistaCCost=        fistaCCost+1/K*showResult(FISTAC012(:,:,k),2,'cost');
        fistaCost=        fistaCost+1/K*showResult(FISTA012(:,:,k),2,'cost');
        spiralCost=        spiralCost+1/K*showResult(spiral012(:,:,k),2,'cost');
        fpcasCost=        fpcasCost+1/K*showResult(fpcas012(:,:,k),2,'f');

        npgRMSE=        npgRMSE+1/K*showResult(npg012(:,:,k),2,'RMSE');
        npgCRMSE=        npgCRMSE+1/K*showResult(npgC012(:,:,k),2,'RMSE');
        fistaCRMSE=        fistaCRMSE+1/K*showResult(FISTAC012(:,:,k),2,'RMSE');
        fistaRMSE=        fistaRMSE+1/K*showResult(FISTA012(:,:,k),2,'RMSE');
        spiralRMSE=        spiralRMSE+1/K*showResult(spiral012(:,:,k),2,'reconerror');
        fpcasRMSE=        fpcasRMSE+1/K*showResult(fpcas012(:,:,k),2,'RMSE');
    end

    % npgTime(:,1:2)=[];
    % npgCTime(:,1:2)=[];
    % fistaCTime(:,1:2)=[];
    % spiralTime(:,1:2)=[];
    % fpcasTime(:,1:2)=[];
    % npgCost(:,1:2)=[];
    % npgCCost(:,1:2)=[];
    % fistaCCost(:,1:2)=[];
    % spiralCost(:,1:2)=[];
    % fpcasCost(:,1:2)=[];
    % npgRMSE(:,1:2)=[];
    % npgCRMSE(:,1:2)=[];
    % fistaCRMSE(:,1:2)=[];
    % spiralRMSE(:,1:2)=[];
    % fpcasRMSE(:,1:2)=[];
    % fpcasRMSE_(:,1:2)=[];
    % fistaCRMSE_(:,1:2)=[];

    [r,c1]=find(npgRMSE  == repmat(min(npgRMSE   ,[],2),1,5)); [r,idx1]=sort(r); [r,c1(idx1)]
    [r,c2]=find(npgCRMSE == repmat(min(npgCRMSE  ,[],2),1,5)); [r,idx2]=sort(r); [r,c2(idx2)]
    [r,c3]=find(fistaCRMSE==repmat(min(fistaCRMSE,[],2),1,5)); [r,idx3]=sort(r); [r,c3(idx3)]
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r); [r,c4(idx4)]
    [r,c5]=find(fpcasRMSE== repmat(min(fpcasRMSE ,[],2),1,5)); [r,idx5]=sort(r); [r,c5(idx5)]
    [r,c6]=find(fistaRMSE== repmat(min(fistaRMSE ,[],2),1,5)); [r,idx6]=sort(r); [r,c6(idx6)]

    figure;
    semilogy(m,npgRMSE   ((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,npgCRMSE  ((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,fistaCRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'r-^');
    semilogy(m,fpcasRMSE ((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m,fistaRMSE ((c6(idx6)-1)*9+(1:9)'),'b-.');
    figure;
    semilogy(m,npgTime   ((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,npgCTime  ((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,fistaCTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'r-^');
    semilogy(m,fpcasTime ((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m,fistaTime ((c6(idx6)-1)*9+(1:9)'),'b-.');
    
    keyboard

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
    forSave=[forSave, sum(fistaCRMSE_,2)./sum(fistaCRMSE_>0,2)];
    forSave=[forSave, sum(fpcasRMSE_,2)./sum(fpcasRMSE_>0,2)];
    save('varyMeasurement.data','forSave','-ascii');
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
    snr=[1 2 5 10 20 50 100 200 500 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10];
    u=[10,1,0.1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];

    for i=(length(snr)-3):length(snr)
        opt.m=m; opt.snr=snr(i);
        opt=loadLinear(conf,opt);
        for j=8:10
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
    snr=[10,  50, 100, 200,  500,  1e3,  1e4,  1e5,  1e6, 1e7];
    u=  [ 1, 0.1, 0.1, 0.1, 1e-2, 1e-2, 1e-2, 1e-3, 1e-4, 1e-5]; % this is for m=600
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
                'reconerrortype',3,...
                'savesolutionpath',0,'verbose',100);
            out.opt=opt; spiral005{i,j}=out;
            spiral005{i,j}.fVal(1)=0.5*norm(conf.Phi(spiral005{i,j}.alpha)-conf.y)^2;
            spiral005{i,j}.fVal(2)=norm(spiral005{i,j}.alpha.*(spiral005{i,j}.alpha<0))^2;
            spiral005{i,j}.fVal(3)=sum(abs(conf.Psit(spiral005{i,j}.alpha)));
            save(filename,'spiral005','-append');

            A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
            option.mxitr=opt.maxItr;
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            fpcas005{i,j}=out; fpcas005{i,j}.alpha = conf.Psi(s);
            fpcas005{i,j}.fVal(1)=0.5*norm(conf.Phi(fpcas005{i,j}.alpha)-conf.y)^2;
            fpcas005{i,j}.fVal(2)=norm(fpcas005{i,j}.alpha.*(fpcas005{i,j}.alpha<0))^2;
            fpcas005{i,j}.fVal(3)=sum(abs(conf.Psit(fpcas005{i,j}.alpha)));
            fpcas005{i,j}.opt = opt; alphaHat=fpcas005{i,j}.alpha;
            fpcas005{i,j}.RMSE=(norm(alphaHat-opt.trueAlpha)/norm(opt.trueAlpha))^2;
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
                'reconerrortype',3,...
                'savesolutionpath',0,'verbose',100);
            out.opt=opt; spiral015{i,j}=out;
            spiral015{i,j}.fVal(1)=0.5*norm(conf.Phi(spiral015{i,j}.alpha)-conf.y)^2;
            spiral015{i,j}.fVal(2)=norm(spiral015{i,j}.alpha.*(spiral015{i,j}.alpha<0))^2;
            spiral015{i,j}.fVal(3)=sum(abs(conf.Psit(spiral015{i,j}.alpha)));
            save(filename,'spiral015','-append');

            A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
            option.mxitr=opt.maxItr;
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            fpcas015{i,j}=out; fpcas015{i,j}.alpha = conf.Psi(s);
            fpcas015{i,j}.fVal(1)=0.5*norm(conf.Phi(fpcas015{i,j}.alpha)-conf.y)^2;
            fpcas015{i,j}.fVal(2)=norm(fpcas015{i,j}.alpha.*(fpcas015{i,j}.alpha<0))^2;
            fpcas015{i,j}.fVal(3)=sum(abs(conf.Psit(fpcas015{i,j}.alpha)));
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
    opt.matrixType='nonneg';
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

        % opt.alphaStep='IST_ADMM_NNL1';
        % ist006{i,j}=lasso(conf.Phi,conf.Phit,...
        %     conf.Psi,conf.Psit,conf.y,initSig,opt);
        % save(filename,'ist006','-append');

        % subtolerance=1e-5;
        % [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
        %     SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
        %     'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','poisson',...
        %     'initialization',initSig,'maxiter',opt.maxItr,...
        %     'miniter',0,'stopcriterion',3,...
        %     'tolerance',opt.thresh,'truth',opt.trueAlpha,...
        %     'subtolerance',subtolerance,'monotone',1,...
        %     'saveobjective',1,'savereconerror',1,'savecputime',1,...
        %     'reconerrortype',3,...
        %     'savesolutionpath',0,'verbose',100);
        % out.opt=opt; spiral006{i,j}=out;
        % save(filename,'spiral006','-append');
    end
end

% The X-ray CT example, test and find the best u for each prjFull
if(any(runList==007))
    load(filename,'*007');
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u=10.^[-1 -2 -3 -4 -5 -6];
    opt.maxItr=2e3; opt.thresh=1e-12;
    for i=1:length(prjFull)
        for j=5:6
            fprintf('%s, i=%d, j=%d\n','X-ray CT example',i,j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2; opt.u = u(j);
            opt=conf.setup(opt);
            conf.y=conf.Phi(opt.trueAlpha); % equivalent to linear projection
            initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

            out007{i}.img=conf.FBP(conf.y);
            out007{i}.alpha=out007{i}.img(opt.mask~=0);
            out007{i}.RSE=norm(conf.y-conf.Phi(out007{i}.alpha))/norm(conf.y);
            out007{i}.RMSE=1-(out007{i}.alpha'*opt.trueAlpha/norm(out007{i}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out007','-append');

            opt.alphaStep='FISTA_ADMM_NNL1';
            out007{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out007','-append');
        end
    end
end

% The X-ray CT example, test and find the best u for each prjFull
if(any(runList==008))     % FPCAS
    load(filename,'*008');
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u=10.^[-6 -5 -5 -5 -5 -5];
    opt.maxItr=2e3; opt.thresh=1e-12;
    for i=1:length(prjFull)
        fprintf('%s, i=%d, j=%d\n','X-ray CT example',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2; opt.u = u(i);
        opt=conf.setup(opt);
        conf.y=conf.Phi(opt.trueAlpha); % equivalent to linear projection
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

        %fbp008{i}.img=conf.FBP(conf.y);
        %fbp008{i}.alpha=fbp008{i}.img(opt.mask~=0);
        %fbp008{i}.RSE=(norm(conf.y-conf.Phi(fbp008{i}.alpha))/norm(conf.y))^2;
        %%fbp008{i}.RMSE=1-(fbp008{i}.alpha'*opt.trueAlpha/norm(fbp008{i}.alpha)/norm(opt.trueAlpha))^2;
        %fbp008{i}.RMSE=(norm(fbp008{i}.alpha-opt.trueAlpha)/norm(opt.trueAlpha))^2;
        %fprintf('fbp RMSE=%g\n',fbp008{i}.RMSE);
        %save(filename,'fbp008','-append');

        if(i>1)
        opt.alphaStep='FISTA_ADMM_NNL1';
        npg008{i,j}=lasso(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'npg008','-append');
        end

        % A = @(xx) conf.Phi(conf.Psi(xx));
        % At = @(yy) conf.Psit(conf.Phit(yy));
        % AO=A_operator(A,At);
        % mu=opt.u; option.x0=conf.Psit(initSig);
        % option.mxitr=opt.maxItr;
        % [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
        % fpcas008{i,j}=out; fpcas008{i,j}.alpha = conf.Psi(s);
        % fpcas008{i,j}.opt = opt;
        % fpcas008{i,j}.fVal(1)=0.5*norm(conf.Phi(fpcas008{i,j}.alpha)-conf.y)^2;
        % fpcas008{i,j}.fVal(2)=norm(fpcas008{i,j}.alpha.*(fpcas008{i,j}.alpha<0))^2;
        % fpcas008{i,j}.fVal(3)=sum(abs(conf.Psit(fpcas008{i,j}.alpha)));
        % fpcas008{i,j}.RSE=(norm(conf.y-conf.Phi(fpcas008{i,j}.alpha))/norm(conf.y))^2;
        % %fpcas008{i,j}.RMSE=1-(fpcas008{i,j}.alpha'*opt.trueAlpha/norm(fpcas008{i,j}.alpha)/norm(opt.trueAlpha))^2;
        % fpcas008{i,j}.RMSE=(norm(fpcas008{i,j}.alpha-opt.trueAlpha)/norm(opt.trueAlpha))^2;
        % fprintf('fpcas RMSE=%g\n',fpcas008{i}.RMSE);
        % save(filename,'fpcas008','-append');

        out=[]; subtolerance=1e-5;
        [out.alpha, out.p, out.cost, out.reconerror, out.time] = ...
            SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',3,...
            'savesolutionpath',0,'verbose',100);
        out.opt=opt; spiral008{i,j}=out;
        spiral008{i,j}.fVal(1)=0.5*norm(conf.Phi(spiral008{i,j}.alpha)-conf.y)^2;
        spiral008{i,j}.fVal(2)=norm(spiral008{i,j}.alpha.*(spiral008{i,j}.alpha<0))^2;
        spiral008{i,j}.fVal(3)=sum(abs(conf.Psit(spiral008{i,j}.alpha)));
        save(filename,'spiral008','-append');
    end
end

% Plot the figures, or save the data for gnuplot in the paper
if(any(runList==999))
    keyboard;

    clear *Time *Cost *RMSE forSave
    snr=[10 50 100 200 500 1e3 1e4 1e5 1e6 1e7];

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
    fistaCRMSE_=showResult(FISTAC005,4,2);
    fpcasRMSE_=showResult(fpcas005,4,2);

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
    forSave=[forSave, sum(fistaCRMSE_,2)./sum(fistaCRMSE_>0,2)];
    forSave=[forSave, sum(fpcasRMSE_,2)./sum(fpcasRMSE_>0,2)];
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
    idx=(forSave(:,1)~=0); forSave(idx,1)=(forSave(idx,1)-mincost);
    idx=(forSave(:,4)~=0); forSave(idx,4)=(forSave(idx,4)-mincost);
    idx=(forSave(:,7)~=0); forSave(idx,7)=(forSave(idx,7)-mincost);
    save('cost_itr.data','forSave','-ascii');

    keyboard
    !cp vary*.data cost_itr.data skyline.data ~/research/myPaper/asilomar2014/
    return;

    i=3; t=0;
    prjFull = [60, 80, 100, 120, 180, 360];
    perform=[]; forSave=[];
    load(filename,'npg008'); out=npg008{i};
    fprintf('NPGwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'NPGwLin.png','png');
    temp=showResult(npg008,2,'RMSE');
    perform=[perform, temp];
    %perform=[perform, [temp(1,6); temp(2:end,5)]];
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    load(filename,'spiral008'); out=spiral008{i};
    fprintf('SPIRALwLin: %g\n',out.reconerror(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'SPIRALwLin.png','png');
    temp=showResult(spiral008,2,'reconerror');
    perform=[perform, temp(:)];
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.reconerror),t)=out.reconerror;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    load('runAsilomar2014','fbp008'); out=fbp008{i};
    fprintf('FBPwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,spiral008{i}.opt.mask);
    imwrite(img/max(img(:)),'FBPwLin.png','png');
    temp=showResult(fbp008,2,'RMSE');
    perform=[perform, temp(:)];

    load('runAsilomar2014','fpcas008');
    out=fpcas008{i};
    fprintf('FPCASwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'FPCASwLin.png','png');
    temp=showResult(fpcas008,2,'RMSE');
    perform=[perform, temp(:,j)];

    perform=[perform,prjFull(:)];
    save('varyPrj.data','perform','-ascii');
    save('xray_itr.data','forSave','-ascii');
    !cp varyPrj.data xray_itr.data *wLin.png ~/research/myPaper/asilomar2014/

    clear *Time *Cost *RMSE forSave
end
end

