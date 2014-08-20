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
paperDir = '~/research/myPaper/asilomar2014/';

% runList rules, abc
% a:
%       0: linear example
%       1: wrist example

%%%%%%%%%%%%%%%%%%%%%%%%

if(any(runList==001))
    filename = [mfilename '_002.mat'];
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
    u = [1e-2,1e-2,1e-3,1e-3,1e-4,1e-4,1e-5,1e-5,1e-5];
    for k=1:20
        for i=1:length(m)
            opt.m=m(i); opt.snr=inf;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            if(k<11) continue; end;
            for j=2:4
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);

                opt.u = u(i)*10^(j-3);

                opt.continuation=false;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npg{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npg','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgc{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgc','-append');

                A = @(xx) conf.Phi(conf.Psi(xx));
                At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; 
                option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
                [s, out] = FPC_AS_mod(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas{i,j,k}=out; fpcas{i,j,k}.alpha = conf.Psi(s);
                fpcas{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(fpcas{i,j,k}.alpha)-conf.y);
                fpcas{i,j,k}.fVal(2)=sqrNorm(fpcas{i,j,k}.alpha.*(fpcas{i,j,k}.alpha<0));
                fpcas{i,j,k}.fVal(3)=pNorm(conf.Psit(fpcas{i,j,k}.alpha),1);
                fpcas{i,j,k}.opt = opt; alphaHat=fpcas{i,j,k}.alpha;
                fpcas{i,j,k}.RMSE=sqrNorm(alphaHat-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
                fprintf('fpcas RMSE=%g\n',fpcas{i,j,k}.RMSE);
                save(filename,'fpcas','-append');

                opt.continuation=false; opt.alphaStep='FISTA_L1';
                npgs{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgs','-append');

                opt.continuation=false; opt.alphaStep='FISTA_L1';
                opt.initStep='fixed'; opt.adaptiveStep=false;
                fista{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'fista','-append');
                opt=rmfield(opt,'initStep'); opt=rmfield(opt,'adaptiveStep');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,'savedifalpha',1,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral{i,j,k}=out;
                spiral{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral{i,j,k}.alpha)-conf.y);
                spiral{i,j,k}.fVal(2)=sqrNorm(spiral{i,j,k}.alpha.*(spiral{i,j,k}.alpha<0));
                spiral{i,j,k}.fVal(3)=pNorm(conf.Psit(spiral{i,j,k}.alpha),1);
                save(filename,'spiral','-append');

                ppsi = @(yyy,uuu) conf.Psi(Utils.softThresh(conf.Psit(yyy),uuu));
                rrrr = @(xxx) pNorm(conf.Psit(xxx),1);
                [x_SpaRSA,x_debias_SpaRSA,obj_SpaRSA,times_SpaRSA,debias_start_SpaRSA,mse]=...
                    SpaRSA_mod(conf.y,conf.Phi,opt.u,...
                    'AT',conf.Phit,...
                    'Psi',ppsi,...
                    'Phi',rrrr,...
                    'Initialization',initSig,...
                    'StopCriterion',5,...
                    'ToleranceA',opt.thresh, ...
                    'True_x',opt.trueAlpha,...
                    'BB_variant',1,...
                    'Safeguard',1,...
                    'Monotone',0,...
                    'Continuation',1,...
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
end

if(any(runList==902))
    filename = [mfilename '_002.mat']; load(filename);

    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    K = 10;

    npgTime   = 0;
    npgcTime  = 0;
    npgsTime  = 0;
    spiralTime= 0;
    fpcasTime = 0;
    fistaTime = 0;
    sparsaTime= 0;
    npgCost   = 0;
    npgcCost  = 0;
    npgsCost  = 0;
    spiralCost= 0;
    fpcasCost = 0;
    fistaCost = 0;
    sparsaCost= 0;
    npgRMSE   = 0;
    npgcRMSE  = 0;
    npgsRMSE  = 0;
    spiralRMSE= 0;
    fpcasRMSE = 0;
    fistaRMSE = 0;
    sparsaRMSE= 0;

    for k=1:K
        npgTime   =   npgTime+1/K*showResult(   npg(:,:,k),2,'time');
        npgcTime  =  npgcTime+1/K*showResult(  npgc(:,:,k),2,'time');
        npgsTime  =  npgsTime+1/K*showResult(  npgs(:,:,k),2,'time');
        spiralTime=spiralTime+1/K*showResult(spiral(:,:,k),2,'time');
        fpcasTime = fpcasTime+1/K*showResult( fpcas(:,:,k),2,'cpu');
        fistaTime = fistaTime+1/K*showResult( fista(:,:,k),2,'time');
        sparsaTime=sparsaTime+1/K*showResult(sparsa(:,:,k),2,'time');

        npgCost   =   npgCost+1/K*showResult(   npg(:,:,k),2,'cost');
        npgcCost  =  npgcCost+1/K*showResult(  npgc(:,:,k),2,'cost');
        npgsCost  =  npgsCost+1/K*showResult(  npgs(:,:,k),2,'cost');
        spiralCost=spiralCost+1/K*showResult(spiral(:,:,k),2,'cost');
        fpcasCost = fpcasCost+1/K*showResult( fpcas(:,:,k),2,'f');
        fistaCost = fistaCost+1/K*showResult( fista(:,:,k),2,'cost');
        sparsaCost=sparsaCost+1/K*showResult(sparsa(:,:,k),2,'cost');

        npgRMSE   =   npgRMSE+1/K*showResult(   npg(:,:,k),2,'RMSE');
        npgcRMSE  =  npgcRMSE+1/K*showResult(  npgc(:,:,k),2,'RMSE');
        npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs(:,:,k),2,'RMSE');
        spiralRMSE=spiralRMSE+1/K*showResult(spiral(:,:,k),2,'reconerror');
        fpcasRMSE = fpcasRMSE+1/K*showResult( fpcas(:,:,k),2,'RMSE');
        fistaRMSE = fistaRMSE+1/K*showResult( fista(:,:,k),2,'RMSE');
        sparsaRMSE=sparsaRMSE+1/K*showResult(sparsa(:,:,k),2,'RMSE');

        % npgRMSE   =   npgRMSE+1/K*showResult(   npg(:,:,k),4,2);
        % npgcRMSE  =  npgcRMSE+1/K*showResult(  npgc(:,:,k),4,2);
        % npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs(:,:,k),4,2);
        % spiralRMSE=spiralRMSE+1/K*showResult(spiral(:,:,k),4,2);
        % fpcasRMSE = fpcasRMSE+1/K*showResult( fpcas(:,:,k),4,2);
        % fistaRMSE = fistaRMSE+1/K*showResult( fista(:,:,k),4,2);
        % sparsaRMSE=sparsaRMSE+1/K*showResult(sparsa(:,:,k),4,2);
    end

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c3]=find(  npgsRMSE== repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c4]=find(spiralRMSE== repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r);
    [r,c7]=find(sparsaRMSE== repmat(min(sparsaRMSE,[],2),1,5)); [r,idx7]=sort(r);
    [c1(idx1) ,c2(idx2) ,c3(idx3) ,c4(idx4) ,c5(idx5) ,c6(idx6) ,c7(idx7)]
    for k=1:10; for i=1:size(npg,1)
        figure(1); semilogy(m(i),   npg{i,gEle(c1(idx1),i),k}.RMSE(end),'r-*');     hold on;
        figure(2);semilogy(m(i),  npgc{i,gEle(c2(idx2),i),k}.RMSE(end),'c-p');hold on;
        figure(3);semilogy(m(i),  npgs{i,gEle(c3(idx3),i),k}.RMSE(end),'k-s');hold on;
        figure(4);semilogy(m(i),spiral{i,gEle(c4(idx4),i),k}.reconerror(end),'k-^');hold on;
        figure(5);semilogy(m(i), fpcas{i,gEle(c5(idx5),i),k}.RMSE(end),'g-o');hold on;
        figure(6);semilogy(m(i), fista{i,gEle(c6(idx6),i),k}.RMSE(end),'b-.');hold on;
        figure(7);semilogy(m(i),sparsa{i,gEle(c7(idx7),i),k}.RMSE(end),'y-p');hold on;
    end; end;
    figure;
    semilogy(m,   npgRMSE((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcRMSE((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m, fpcasRMSE((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaRMSE((c6(idx6)-1)*9+(1:9)'),'b-.');
    semilogy(m,sparsaRMSE((c7(idx7)-1)*9+(1:9)'),'y-p');
    figure;
    semilogy(m,   npgTime((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcTime((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m, fpcasTime((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaTime((c6(idx6)-1)*9+(1:9)'),'b-.');
    semilogy(m,sparsaTime((c7(idx7)-1)*9+(1:9)'),'y-p');

    keyboard

    idx1 = 2;    
    idx2 = 2;
    idx3 = 2;
    idx4 = 3;
    idx5 = 2;
    idx6 = 2;
    idx7 = 2;
    figure;
    semilogy(m,   npgRMSE(:,idx1),'r-*'); hold on;
    semilogy(m,  npgcRMSE(:,idx2),'c-p');
    semilogy(m,  npgsRMSE(:,idx3),'k-s');
    semilogy(m,spiralRMSE(:,idx4),'k-^');
    semilogy(m, fpcasRMSE(:,idx5),'g-o');
    semilogy(m, fistaRMSE(:,idx6),'b-.');
    semilogy(m,sparsaRMSE(:,idx7),'y-p');
    figure;
    plot(m,-10*log10(   npgRMSE(:,idx1)),'r-*'); hold on;
    plot(m,-10*log10(  npgcRMSE(:,idx2)),'c-p');
    plot(m,-10*log10(  npgsRMSE(:,idx3)),'k-s');
    plot(m,-10*log10(spiralRMSE(:,idx4)),'k-^');
    plot(m,-10*log10( fpcasRMSE(:,idx5)),'g-o');
    plot(m,-10*log10( fistaRMSE(:,idx6)),'b-.');
    plot(m,-10*log10(sparsaRMSE(:,idx7)),'y-p');
    figure;
    semilogy(m,   npgTime(:,idx1),'r-*'); hold on;
    semilogy(m,  npgcTime(:,idx2),'c-p');
    semilogy(m,  npgsTime(:,idx3),'k-s');
    semilogy(m,spiralTime(:,idx4),'k-^');
    semilogy(m, fpcasTime(:,idx5),'g-o');
    semilogy(m, fistaTime(:,idx6),'b-.');
    semilogy(m,sparsaTime(:,idx7),'y-p');

    keyboard
    
    temp = 4;
    signal=npg{1}.opt.trueAlpha;
    signal=[signal,    npg{temp,idx1,1}.alpha];
    signal=[signal,   npgc{temp,idx2,1}.alpha];
    signal=[signal,   npgs{temp,idx3,1}.alpha];
    signal=[signal, spiral{temp,idx4,1}.alpha];
    signal=[signal,  fpcas{temp,idx5,1}.alpha];
    signal=[signal,  fista{temp,idx6,1}.alpha];
    signal=[signal, sparsa{temp,idx7,1}.alpha];
    save('skyline.data','signal','-ascii');
    figure; plot(signal(:,2)); hold on; plot(signal(:,1),'r'); title('NPG');
    figure; plot(signal(:,4)); hold on; plot(signal(:,1),'r'); title('NPGs');
    figure; plot(signal(:,6)); hold on; plot(signal(:,1),'r'); title('FPCas');

    M=length(m);
    str=        '$m$            ';              for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%10d',str,m(i)); end; end;
    str=sprintf('%s\\\\\\hline',str);
    str=sprintf('%s\nFISTA          ', str);    for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str, fista{i,idx6,1}.cost(end));end; end;
    str=sprintf('%s\\\\\nNPG$_\\text{S}$ ',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,  npgs{i,idx3,1}.cost(end));end; end;
    str=sprintf('%s\\\\\nFPC$_\\text{AS}$',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str, fpcas{i,idx5,1}.f   (end));end; end;
    str=sprintf('%s\\\\\nSPIRAL$_{-4}$  ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,spiral{i,idx4,1}.cost(end));end; end;
    str=sprintf('%s\\\\\nSPIRAL$_{-5}$  ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,spiral{i,idx6,1}.cost(end));end; end;
    str=sprintf('%s\\\\\nNPG            ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,   npg{i,idx1,1}.cost(end));end; end;
  % str=sprintf('%s\\\\\nNPG$_\\text{c}$ ',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,  npgcCost(i,idx2));end; end;
    file=fopen('varyMeasurementTable.tex','w'); fprintf(file,'%s',str); fclose(file);

    % figure;
    % for i=1:M;
    %     semilogy(npgs{i,idx3,1}.stepSize); hold on; semilogy(fista{i,idx6,1}.stepSize,'r:');
    %     semilogy([1,length(fista{i,idx6,1}.RMSE)],ones(1,2)*1/fista{i,idx6,1}.opt.L,'k-.');
    %     hold off;
    %     pause;
    % end

    forSave=[];
    forSave=[forSave,    npgTime(:,idx1)];
    forSave=[forSave,   npgcTime(:,idx2)];
    forSave=[forSave,   npgsTime(:,idx3)];
    forSave=[forSave, spiralTime(:,idx4)];
    forSave=[forSave,  fpcasTime(:,idx5)];
    forSave=[forSave,  fistaTime(:,idx6)];

    forSave=[forSave,    npgCost(:,idx1)];
    forSave=[forSave,   npgcCost(:,idx2)];
    forSave=[forSave,   npgsCost(:,idx3)];
    forSave=[forSave, spiralCost(:,idx4)];
    forSave=[forSave,  fpcasCost(:,idx5)];
    forSave=[forSave,  fistaCost(:,idx6)];

    forSave=[forSave,    npgRMSE(:,idx1)];
    forSave=[forSave,   npgcRMSE(:,idx2)];
    forSave=[forSave,   npgsRMSE(:,idx3)];
    forSave=[forSave, spiralRMSE(:,idx4)];
    forSave=[forSave,  fpcasRMSE(:,idx5)];
    forSave=[forSave,  fistaRMSE(:,idx6)];
    forSave=[forSave, m(:)];
    forSave=[forSave, spiralTime(:,idx6)];
    forSave=[forSave, spiralCost(:,idx6)];
    forSave=[forSave, spiralRMSE(:,idx6)];
    forSave=[forSave, sparsaTime(:,idx7)];
    forSave=[forSave, sparsaCost(:,idx7)];
    forSave=[forSave, sparsaRMSE(:,idx7)];
    save('varyMeasurement.data','forSave','-ascii');

    mIdx=6; forSave=[]; t=0;
    t=t+1; temp=  npgs{mIdx,idx3,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= fista{mIdx,idx6,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgs{mIdx,idx3,1}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fista{mIdx,idx6,1}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgs{mIdx,idx3,1}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fista{mIdx,idx6,1}.time(:);     forSave(1:length(temp),t)=temp;
    save('stepSizeLin.data','forSave','-ascii');

    system(['mv stepSizeLin.data varyMeasurement.data skyline.data varyMeasurementTable.tex ' paperDir]);
    disp('done');
end

% vary the number of measurements, with continuation, nonneg Phi
if(any(runList==012))
    load(filename,'*012');
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT(); opt.debugLevel=0;
    opt.maxItr=1e4; opt.thresh=1e-6; opt.matrixType='gaussian';
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    for k=1:1
        for i=1:length(m)
            opt.m=m(i); opt.snr=1e6;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = (1e-5)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=1:5
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
                opt.u = u*10^(j-2);

                % opt.continuation=true;
                % opt.alphaStep='FISTA_ADMM_NNL1';
                % npgc012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'npgc012','-append');

                % opt.continuation=false;
                % opt.alphaStep='FISTA_ADMM_NNL1';
                % npg012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'npg012','-append');

                % opt.continuation=true;
                % opt.alphaStep='FISTA_L1';
                % fistac012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'fistac012','-append');

                % opt.continuation=false;
                % opt.alphaStep='FISTA_L1';
                % fista012{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'fista012','-append');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,'savedifalpha',1,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral012{i,j,k}=out;
                spiral012{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral012{i,j,k}.alpha)-conf.y);
                spiral012{i,j,k}.fVal(2)=sqrNorm(spiral012{i,j,k}.alpha.*(spiral012{i,j,k}.alpha<0));
                spiral012{i,j,k}.fVal(3)=sum(abs(conf.Psit(spiral012{i,j,k}.alpha)));
                save(filename,'spiral012','-append');

                A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
                [s, out] = FPC_AS_mod(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas012{i,j,k}=out; fpcas012{i,j,k}.alpha = conf.Psi(s);
                fpcas012{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi( fpcas012{i,j,k}.alpha)-conf.y);
                fpcas012{i,j,k}.fVal(2)=sqrNorm(fpcas012{i,j,k}.alpha.*(fpcas012{i,j,k}.alpha<0));
                fpcas012{i,j,k}.fVal(3)=sum(abs(conf.Psit(fpcas012{i,j,k}.alpha)));
                fpcas012{i,j,k}.opt = opt; alphaHat=fpcas012{i,j,k}.alpha;
                fpcas012{i,j,k}.RMSE=sqrNorm(alphaHat-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
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

    npgTime   = 0;
    npgcTime  = 0;
    npgsTime  = 0;
    fistaTime = 0;
    spiralTime= 0;
    fpcasTime = 0;
    npgCost   = 0;
    npgcCost  = 0;
    npgsCost  = 0;
    fistaCost = 0;
    spiralCost= 0;
    fpcasCost = 0;
    npgRMSE   = 0;
    npgcRMSE  = 0;
    npgsRMSE  = 0;
    fistaRMSE = 0;
    spiralRMSE= 0;
    fpcasRMSE = 0;

    for k=1:K
        npgTime   =   npgTime+1/K*showResult(   npg012(:,:,k),2,'time');
        npgcTime  =  npgcTime+1/K*showResult(  npgc012(:,:,k),2,'time');
        npgsTime  =  npgsTime+1/K*showResult(  npgs012(:,:,k),2,'time');
        fistaTime = fistaTime+1/K*showResult( fista012(:,:,k),2,'time');
        spiralTime=spiralTime+1/K*showResult(spiral012(:,:,k),2,'time');
        fpcasTime = fpcasTime+1/K*showResult( fpcas012(:,:,k),2,'cpu');

        npgCost   =   npgCost+1/K*showResult(   npg012(:,:,k),2,'cost');
        npgcCost  =  npgcCost+1/K*showResult(  npgc012(:,:,k),2,'cost');
        npgsCost  =  npgsCost+1/K*showResult(  npgs012(:,:,k),2,'cost');
        fistaCost = fistaCost+1/K*showResult( fista012(:,:,k),2,'cost');
        spiralCost=spiralCost+1/K*showResult(spiral012(:,:,k),2,'cost');
        fpcasCost = fpcasCost+1/K*showResult( fpcas012(:,:,k),2,'f');

        npgRMSE   =   npgRMSE+1/K*showResult(   npg012(:,:,k),2,'RMSE');
        npgcRMSE  =  npgcRMSE+1/K*showResult(  npgc012(:,:,k),2,'RMSE');
        npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs012(:,:,k),2,'RMSE');
        fistaRMSE = fistaRMSE+1/K*showResult( fista012(:,:,k),2,'RMSE');
        spiralRMSE=spiralRMSE+1/K*showResult(spiral012(:,:,k),2,'reconerror');
        fpcasRMSE = fpcasRMSE+1/K*showResult( fpcas012(:,:,k),2,'RMSE');
    end

    % npgTime(:,1:2)=[];
    % npgcTime(:,1:2)=[];
    %  npgsTime(:,1:2)=[];
    % spiralTime(:,1:2)=[];
    % fpcasTime(:,1:2)=[];
    % npgCost(:,1:2)=[];
    % npgcCost(:,1:2)=[];
    %  npgsCost(:,1:2)=[];
    % spiralCost(:,1:2)=[];
    % fpcasCost(:,1:2)=[];
    % npgRMSE(:,1:2)=[];
    % npgcRMSE(:,1:2)=[];
    %  npgsRMSE(:,1:2)=[];
    % spiralRMSE(:,1:2)=[];
    % fpcasRMSE(:,1:2)=[];
    % fpcasRMSE_(:,1:2)=[];
    %  npgsRMSE_(:,1:2)=[];

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r); [r,c1(idx1)]
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r); [r,c2(idx2)]
    [r,c3]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r); [r,c3(idx3)]
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r); [r,c4(idx4)]
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r); [r,c5(idx5)]
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r); [r,c6(idx6)]

    figure;
    semilogy(m,   npgRMSE((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcRMSE((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'r-^');
    semilogy(m, fpcasRMSE((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaRMSE((c6(idx6)-1)*9+(1:9)'),'b-.');
    figure;
    semilogy(m,   npgTime((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcTime((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'r-^');
    semilogy(m, fpcasTime((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaTime((c6(idx6)-1)*9+(1:9)'),'b-.');

    figure;
    semilogy(m,   npgRMSE(:,3),'r-*'); hold on;
    semilogy(m,  npgcRMSE(:,3),'c-p');
    semilogy(m,  npgsRMSE(:,3),'k-s');
    semilogy(m,spiralRMSE(:,3),'k-^');
    semilogy(m, fpcasRMSE(:,2),'g-o');
    semilogy(m, fpcasRMSE(:,3),'g-<');
    semilogy(m, fistaRMSE(:,3),'b-.');
    figure;
    semilogy(m,   npgTime(:,3),'r-*'); hold on;
    semilogy(m,  npgcTime(:,3),'c-p');
    semilogy(m,  npgsTime(:,3),'k-s');
    semilogy(m,spiralTime(:,3),'k-^');
    semilogy(m, fpcasTime(:,2),'g-o');
    semilogy(m, fpcasTime(:,3),'g-<');
    semilogy(m, fistaTime(:,3),'b-.');


    keyboard

    forSave=[];
    forSave=[forSave, sum(   npgTime,2)./sum(   npgTime>0,2)];
    forSave=[forSave, sum(  npgcTime,2)./sum(  npgcTime>0,2)];
    forSave=[forSave, sum(  npgsTime,2)./sum(  npgsTime>0,2)];
    forSave=[forSave, sum(spiralTime,2)./sum(spiralTime>0,2)];
    forSave=[forSave, sum( fpcasTime,2)./sum( fpcasTime>0,2)];

    forSave=[forSave, sum(   npgCost,2)./sum(   npgCost>0,2)];
    forSave=[forSave, sum(  npgcCost,2)./sum(  npgcCost>0,2)];
    forSave=[forSave, sum(  npgsCost,2)./sum(  npgsCost>0,2)];
    forSave=[forSave, sum(spiralCost,2)./sum(spiralCost>0,2)];
    forSave=[forSave, sum( fpcasCost,2)./sum( fpcasCost>0,2)];

    forSave=[forSave, sum(   npgRMSE,2)./sum(   npgRMSE>0,2)];
    forSave=[forSave, sum(  npgcRMSE,2)./sum(  npgcRMSE>0,2)];
    forSave=[forSave, sum(  npgsRMSE,2)./sum(  npgsRMSE>0,2)];
    forSave=[forSave, sum(spiralRMSE,2)./sum(spiralRMSE>0,2)];
    forSave=[forSave, sum( fpcasRMSE,2)./sum( fpcasRMSE>0,2)];
    forSave=[forSave, m(:)];
    forSave=[forSave, sum(  npgsRMSE_,2)./sum(  npgsRMSE_>0,2)];
    forSave=[forSave, sum( fpcasRMSE_,2)./sum( fpcasRMSE_>0,2)];
    save('varyMeasurement.data','forSave','-ascii');
end

% vary the number of measurements, with continuation, nonneg Phi
if(any(runList==022))
    load(filename,'*022');
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT(); opt.debugLevel=0;
    opt.maxItr=1e4; opt.thresh=1e-6; opt.matrixType='gaussian';
    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    for k=1:1
        for i=1:length(m)
            opt.m=m(i); opt.snr=1e7;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = (1e-5)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=1:5
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);
                opt.u = u*10^(j-2);

                % opt.continuation=true;
                % opt.alphaStep='FISTA_ADMM_NNL1';
                % npgc022{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'npgc022','-append');

                % opt.continuation=false;
                % opt.alphaStep='FISTA_ADMM_NNL1';
                % npg022{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'npg022','-append');

                % opt.continuation=true;
                % opt.alphaStep='FISTA_L1';
                % fistac022{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'fistac022','-append');

                % opt.continuation=false;
                % opt.alphaStep='FISTA_L1';
                % fista022{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'fista022','-append');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,'savedifalpha',1,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral022{i,j,k}=out;
                spiral022{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral022{i,j,k}.alpha)-conf.y);
                spiral022{i,j,k}.fVal(2)=sqrNorm(spiral022{i,j,k}.alpha.*(spiral022{i,j,k}.alpha<0));
                spiral022{i,j,k}.fVal(3)=sum(abs(conf.Psit(spiral022{i,j,k}.alpha)));
                save(filename,'spiral022','-append');

                A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
                [s, out] = FPC_AS_mod(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas022{i,j,k}=out; fpcas022{i,j,k}.alpha = conf.Psi(s);
                fpcas022{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi( fpcas022{i,j,k}.alpha)-conf.y);
                fpcas022{i,j,k}.fVal(2)=sqrNorm(fpcas022{i,j,k}.alpha.*(fpcas022{i,j,k}.alpha<0));
                fpcas022{i,j,k}.fVal(3)=sum(abs(conf.Psit(fpcas022{i,j,k}.alpha)));
                fpcas022{i,j,k}.opt = opt; alphaHat=fpcas022{i,j,k}.alpha;
                fpcas022{i,j,k}.RMSE=sqrNorm(alphaHat-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
                fprintf('fpcas RMSE=%g\n',fpcas022{i,j,k}.RMSE);
                save(filename,'fpcas022','-append');
            end
        end
    end
end

if(any(runList==922))
    load(filename,'*022'); j=1;
    opt.a=[];
    opt=loadLinear(ConfigCT(),opt);
    signal=opt.trueAlpha;
    save('skyline.data','signal','-ascii');

    m=[ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    K = 1;

npgTime   = 0;
npgcTime  = 0;
npgsTime  = 0;
fistaTime = 0;
spiralTime= 0;
fpcasTime = 0;
npgCost   = 0;
npgcCost  = 0;
npgsCost  = 0;
fistaCost = 0;
spiralCost= 0;
fpcasCost = 0;
npgRMSE   = 0;
npgcRMSE  = 0;
npgsRMSE  = 0;
fistaRMSE = 0;
spiralRMSE= 0;
fpcasRMSE = 0;

    for k=1:K
npgTime   =   npgTime+1/K*showResult(   npg022(:,:,k),2,'time');
npgcTime  =  npgcTime+1/K*showResult(  npgc022(:,:,k),2,'time');
npgsTime  =  npgsTime+1/K*showResult(  npgs022(:,:,k),2,'time');
fistaTime = fistaTime+1/K*showResult( fista022(:,:,k),2,'time');
spiralTime=spiralTime+1/K*showResult(spiral022(:,:,k),2,'time');
fpcasTime = fpcasTime+1/K*showResult( fpcas022(:,:,k),2,'cpu');

npgCost   =   npgCost+1/K*showResult(   npg022(:,:,k),2,'cost');
npgcCost  =  npgcCost+1/K*showResult(  npgc022(:,:,k),2,'cost');
npgsCost  =  npgsCost+1/K*showResult(  npgs022(:,:,k),2,'cost');
fistaCost = fistaCost+1/K*showResult( fista022(:,:,k),2,'cost');
spiralCost=spiralCost+1/K*showResult(spiral022(:,:,k),2,'cost');
fpcasCost = fpcasCost+1/K*showResult( fpcas022(:,:,k),2,'f');

npgRMSE   =   npgRMSE+1/K*showResult(   npg022(:,:,k),2,'RMSE');
npgcRMSE  =  npgcRMSE+1/K*showResult(  npgc022(:,:,k),2,'RMSE');
npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs022(:,:,k),2,'RMSE');
fistaRMSE = fistaRMSE+1/K*showResult( fista022(:,:,k),2,'RMSE');
spiralRMSE=spiralRMSE+1/K*showResult(spiral022(:,:,k),2,'reconerror');
fpcasRMSE = fpcasRMSE+1/K*showResult( fpcas022(:,:,k),2,'RMSE');
    end

    % npgTime(:,1:2)=[];
    % npgcTime(:,1:2)=[];
    %  npgsTime(:,1:2)=[];
    % spiralTime(:,1:2)=[];
    % fpcasTime(:,1:2)=[];
    % npgCost(:,1:2)=[];
    % npgcCost(:,1:2)=[];
    %  npgsCost(:,1:2)=[];
    % spiralCost(:,1:2)=[];
    % fpcasCost(:,1:2)=[];
    % npgRMSE(:,1:2)=[];
    % npgcRMSE(:,1:2)=[];
    %  npgsRMSE(:,1:2)=[];
    % spiralRMSE(:,1:2)=[];
    % fpcasRMSE(:,1:2)=[];
    % fpcasRMSE_(:,1:2)=[];
    %  npgsRMSE_(:,1:2)=[];

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r); [r,c1(idx1)]
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r); [r,c2(idx2)]
    [r,c3]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r); [r,c3(idx3)]
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r); [r,c4(idx4)]
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r); [r,c5(idx5)]
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r); [r,c6(idx6)]

    figure;
    semilogy(m,   npgRMSE((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcRMSE((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'r-^');
    semilogy(m, fpcasRMSE((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaRMSE((c6(idx6)-1)*9+(1:9)'),'b-.');
    figure;
    semilogy(m,   npgTime((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcTime((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'r-^');
    semilogy(m, fpcasTime((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaTime((c6(idx6)-1)*9+(1:9)'),'b-.');

    figure;
    semilogy(m,   npgRMSE(:,2),'r-*'); hold on;
    semilogy(m,  npgcRMSE(:,2),'c-p');
    semilogy(m,  npgsRMSE(:,2),'k-s');
    semilogy(m,spiralRMSE(:,3),'k-^');
    semilogy(m, fpcasRMSE(:,2),'g-o');
    semilogy(m, fpcasRMSE(:,3),'g-<');
    semilogy(m, fistaRMSE(:,2),'b-.');
    
    
    keyboard

    forSave=[];
    forSave=[forSave, sum(   npgTime,2)./sum(   npgTime>0,2)];
    forSave=[forSave, sum(  npgcTime,2)./sum(  npgcTime>0,2)];
    forSave=[forSave, sum(  npgsTime,2)./sum(  npgsTime>0,2)];
    forSave=[forSave, sum(spiralTime,2)./sum(spiralTime>0,2)];
    forSave=[forSave, sum( fpcasTime,2)./sum( fpcasTime>0,2)];

    forSave=[forSave, sum(   npgCost,2)./sum(   npgCost>0,2)];
    forSave=[forSave, sum(  npgcCost,2)./sum(  npgcCost>0,2)];
    forSave=[forSave, sum(  npgsCost,2)./sum(  npgsCost>0,2)];
    forSave=[forSave, sum(spiralCost,2)./sum(spiralCost>0,2)];
    forSave=[forSave, sum( fpcasCost,2)./sum( fpcasCost>0,2)];

    forSave=[forSave, sum(   npgRMSE,2)./sum(   npgRMSE>0,2)];
    forSave=[forSave, sum(  npgcRMSE,2)./sum(  npgcRMSE>0,2)];
    forSave=[forSave, sum(  npgsRMSE,2)./sum(  npgsRMSE>0,2)];
    forSave=[forSave, sum(spiralRMSE,2)./sum(spiralRMSE>0,2)];
    forSave=[forSave, sum( fpcasRMSE,2)./sum( fpcasRMSE>0,2)];
    forSave=[forSave, m(:)];
    forSave=[forSave, sum(  npgsRMSE_,2)./sum(  npgsRMSE_>0,2)];
    forSave=[forSave, sum( fpcasRMSE_,2)./sum( fpcasRMSE_>0,2)];
    save('varyMeasurement.data','forSave','-ascii');
end


% vary the SNR of measurements, with continuation (continuation is good) to 
% find their corresponding good u, m=600;
if(any(runList==004))
    load(filename,'*004');
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-6;
    m=[600];
    snr=[1 2 5 10 20 50 100 200 500 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10];
    u=[10,1,0.1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8];

    for k=1:1
        for i=(length(snr)-3):length(snr)
            opt.m=m; opt.snr=snr(i);
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = (1e-5)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=8:10
                fprintf('%s, i=%d, j=%d\n','FISTA_ADMM_NNL1',i,j);
                opt.u = u*10^(j-2);
                opt.u = u(j);

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgc004{i,j}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgc004','-append');
            end
        end
    end
end

% vary the SNR of measurements, with continuation (continuation is good) to 
% find their corresponding good u, m=700;
if(any(runList==014))
    load(filename,'*014');
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.maxItr=1e4; opt.thresh=1e-6;
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
            npgc014{i,j}=lasso(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgc014','-append');
        end
    end
end

% vary the SNR for m=600, u is picked to give the best result
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

    for k=1:5
        for i=1:length(snr)
            opt.m=m; opt.snr=snr(i);
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = a(i)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=1:4
                opt.u = u*10^(j-3);
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);

                opt.continuation=false;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npg{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npg','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgc{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgc','-append');

                opt.continuation=false; opt.alphaStep='FISTA_L1';
                opt.initStep='fixed'; opt.adaptiveStep=false;
                fista{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'fista','-append');
                opt=rmfield(opt,'initStep'); opt=rmfield(opt,'adaptiveStep');

                opt.continuation=false; opt.alphaStep='FISTA_L1';
                npgs{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgs','-append');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,'savedifalpha',1,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral{i,j,k}=out;
                spiral{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral{i,j,k}.alpha)-conf.y);
                spiral{i,j,k}.fVal(2)=sqrNorm(spiral{i,j,k}.alpha.*(spiral{i,j,k}.alpha<0));
                spiral{i,j,k}.fVal(3)=pNorm(conf.Psit(spiral{i,j,k}.alpha),1);
                save(filename,'spiral','-append');

                A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
                [s, out] = FPC_AS_mod(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas{i,j,k}=out; fpcas{i,j,k}.alpha = conf.Psi(s);
                fpcas{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi( fpcas{i,j,k}.alpha)-conf.y);
                fpcas{i,j,k}.fVal(2)=sqrNorm( fpcas{i,j,k}.alpha.*(fpcas{i,j,k}.alpha<0));
                fpcas{i,j,k}.fVal(3)=pNorm(conf.Psit(fpcas{i,j,k}.alpha),1);
                fpcas{i,j,k}.opt = opt; alphaHat=fpcas{i,j,k}.alpha;
                fpcas{i,j,k}.RMSE=sqrNorm(alphaHat-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
                fprintf('fpcas RMSE=%g\n',fpcas{i,j,k}.RMSE);
                save(filename,'fpcas','-append');
            end
        end
    end
end

if(any(runList==905))
    filename = [mfilename '_005.mat']; load(filename);
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7]';
    K = 5;

    npgTime   = 0;
    npgcTime  = 0;
    npgsTime  = 0;
    fistaTime = 0;
    spiralTime= 0;
    fpcasTime = 0;
    npgCost   = 0;
    npgcCost  = 0;
    npgsCost  = 0;
    fistaCost = 0;
    spiralCost= 0;
    fpcasCost = 0;
    npgRMSE   = 0;
    npgcRMSE  = 0;
    npgsRMSE  = 0;
    fistaRMSE = 0;
    spiralRMSE= 0;
    fpcasRMSE = 0;

    for k=1:K
        npgTime   =   npgTime+1/K*showResult(   npg(:,:,k),2,'time');
        npgcTime  =  npgcTime+1/K*showResult(  npgc(:,:,k),2,'time');
        npgsTime  =  npgsTime+1/K*showResult(  npgs(:,:,k),2,'time');
        spiralTime=spiralTime+1/K*showResult(spiral(:,:,k),2,'time');
        fpcasTime = fpcasTime+1/K*showResult( fpcas(:,:,k),2,'cpu');
        fistaTime = fistaTime+1/K*showResult( fista(:,:,k),2,'time');

        npgCost   =   npgCost+1/K*showResult(   npg(:,:,k),2,'cost');
        npgcCost  =  npgcCost+1/K*showResult(  npgc(:,:,k),2,'cost');
        npgsCost  =  npgsCost+1/K*showResult(  npgs(:,:,k),2,'cost');
        spiralCost=spiralCost+1/K*showResult(spiral(:,:,k),2,'cost');
        fpcasCost = fpcasCost+1/K*showResult( fpcas(:,:,k),2,'f');
        fistaCost = fistaCost+1/K*showResult( fista(:,:,k),2,'cost');

        npgRMSE   =   npgRMSE+1/K*showResult(   npg(:,:,k),2,'RMSE');
        npgcRMSE  =  npgcRMSE+1/K*showResult(  npgc(:,:,k),2,'RMSE');
        npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs(:,:,k),2,'RMSE');
        spiralRMSE=spiralRMSE+1/K*showResult(spiral(:,:,k),2,'reconerror');
        fpcasRMSE = fpcasRMSE+1/K*showResult( fpcas(:,:,k),2,'RMSE');
        fistaRMSE = fistaRMSE+1/K*showResult( fista(:,:,k),2,'RMSE');
    end

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c3]=find(  npgsRMSE== repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c4]=find(spiralRMSE== repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r);
    [c1(idx1) ,c2(idx2) ,c3(idx3) ,c4(idx4) ,c5(idx5) ,c6(idx6)]
    idx1 = 3;    
    idx2 = 3;
    idx3 = 3;
    idx4 = 3;
    idx5 = 3;
    idx6 = 3;

    figure;
    loglog(snr,   npgRMSE(:,idx1),'r-*'); hold on;
    loglog(snr,  npgcRMSE(:,idx2),'c-p');
    loglog(snr,  npgsRMSE(:,idx3),'k-s');
    loglog(snr,spiralRMSE(:,idx4),'k-^');
    loglog(snr, fpcasRMSE(:,idx5),'g-<');
    loglog(snr, fistaRMSE(:,idx6),'b-.');

    figure;
    semilogx(snr,-10*log10(   npgRMSE(:,idx1)),'r-*'); hold on;
    semilogx(snr,-10*log10(  npgcRMSE(:,idx2)),'c-p');
    semilogx(snr,-10*log10(  npgsRMSE(:,idx3)),'k-s');
    semilogx(snr,-10*log10(spiralRMSE(:,idx4)),'k-^');
    semilogx(snr,-10*log10( fpcasRMSE(:,idx5)),'g-o');
    semilogx(snr,-10*log10( fistaRMSE(:,idx6)),'b-.');

    figure;
    loglog(snr,   npgTime(:,idx1),'r-*'); hold on;
    loglog(snr,  npgcTime(:,idx2),'c-p');
    loglog(snr,  npgsTime(:,idx3),'k-s');
    loglog(snr,spiralTime(:,idx4),'k-^');
    loglog(snr, fpcasTime(:,idx5),'g-<');
    loglog(snr, fistaTime(:,idx6),'b-.');

    keyboard

    M=length(snr);
    str=        'SNR            ';              for i=1:M; str=sprintf('%s&%10d',str,snr(i)); end;
    str=sprintf('%s\\\\\\hline',str);
    str=sprintf('%s\nNPG            ', str);    for i=1:M; str=sprintf('%s&%-10.6g',str,   npg{i,idx1,1}.cost(end));end;
  % str=sprintf('%s\\\\\nNPG$_\\text{c}$ ',str);for i=1:M; str=sprintf('%s&%-10.6g',str,  npgc{i,idx2,1}.cost(end));end;
    str=sprintf('%s\\\\\nNPG$_\\text{S}$ ',str);for i=1:M; str=sprintf('%s&%-10.6g',str,  npgs{i,idx3,1}.cost(end));end;
    str=sprintf('%s\\\\\nSPIRAL         ', str);for i=1:M; str=sprintf('%s&%-10.6g',str,spiral{i,idx4,1}.cost(end));end;
    str=sprintf('%s\\\\\nFPC$_\\text{AS}$',str);for i=1:M; str=sprintf('%s&%-10.6g',str, fpcas{i,idx5,1}.f   (end));end;
    str=sprintf('%s\\\\\nFISTA          ', str);for i=1:M; str=sprintf('%s&%-10.6g',str, fista{i,idx6,1}.cost(end));end;
    file=fopen('varySNRTable.tex','w'); fprintf(file,'%s',str); fclose(file);

    % figure;
    % for i=1:M;
    %     semilogy(npgs002{i,idx3,1}.stepSize); hold on; semilogy(fista002{i,idx6,1}.stepSize,'r:');
    %     semilogy([1,length(fista002{i,idx6,1}.RMSE)],ones(1,2)*1/fista002{i,idx6,1}.opt.L,'k-.');
    %     hold off;
    %     pause;
    % end

    forSave=[];
    forSave=[forSave,      npgTime(3:end,idx1)];
    forSave=[forSave,     npgcTime(3:end,idx2)];
    forSave=[forSave,     npgsTime(3:end,idx3)];
    forSave=[forSave,   spiralTime(3:end,idx4)];
    forSave=[forSave,    fpcasTime(3:end,idx5)];
    forSave=[forSave,    fistaTime(3:end,idx6)];

    forSave=[forSave,      npgCost(3:end,idx1)];
    forSave=[forSave,     npgcCost(3:end,idx2)];
    forSave=[forSave,     npgsCost(3:end,idx3)];
    forSave=[forSave,   spiralCost(3:end,idx4)];
    forSave=[forSave,    fpcasCost(3:end,idx5)];
    forSave=[forSave,    fistaCost(3:end,idx6)];

    forSave=[forSave,      npgRMSE(3:end,idx1)];
    forSave=[forSave,     npgcRMSE(3:end,idx2)];
    forSave=[forSave,     npgsRMSE(3:end,idx3)];
    forSave=[forSave,   spiralRMSE(3:end,idx4)];
    forSave=[forSave,    fpcasRMSE(3:end,idx5)];
    forSave=[forSave,    fistaRMSE(3:end,idx6)];
    forSave=[forSave, 10*log10(snr(3:end))];
    save('varySNR.data','forSave','-ascii');
    system(['mv varySNRTable.tex varySNR.data ' paperDir]);
end

% vary the SNR for m=400, u is picked to give the best result
if(any(runList==015))
    filename = [mfilename '_015.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    conf=ConfigCT();
    opt.debugLevel=0;
    opt.maxItr=1e4; opt.thresh=1e-6;
    m=[400];
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7];
    a  =[1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3,1e-4,1e-4,1e-5];
    for k=1:3
        for i=1:length(snr)
            opt.m=m; opt.snr=snr(i);
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0;
            u = a(i)*pNorm(conf.Psit(conf.Phit(conf.y)),inf);
            for j=1:5
                opt.u = u*10^(j-3);
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);

                opt.continuation=false;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npg{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npg','-append');

                opt.continuation=true;
                opt.alphaStep='FISTA_ADMM_NNL1';
                npgc{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgc','-append');

                opt.continuation=false; opt.alphaStep='FISTA_L1';
                opt.initStep='fixed'; opt.adaptiveStep=false;
                fista{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'fista','-append');
                opt=rmfield(opt,'initStep'); opt=rmfield(opt,'adaptiveStep');

                opt.continuation=false; opt.alphaStep='FISTA_L1';
                npgs{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgs','-append');

                subtolerance=1e-5;
                [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                    SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                    'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                    'initialization',initSig,'maxiter',opt.maxItr,...
                    'miniter',0,'stopcriterion',3,...
                    'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                    'subtolerance',subtolerance,'monotone',1,...
                    'saveobjective',1,'savereconerror',1,'savecputime',1,...
                    'reconerrortype',3,'savedifalpha',1,...
                    'savesolutionpath',0,'verbose',1000);
                out.opt=opt; spiral{i,j,k}=out;
                spiral{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral{i,j,k}.alpha)-conf.y);
                spiral{i,j,k}.fVal(2)=sqrNorm(spiral{i,j,k}.alpha.*(spiral{i,j,k}.alpha<0));
                spiral{i,j,k}.fVal(3)=pNorm(conf.Psit(spiral{i,j,k}.alpha),1);
                save(filename,'spiral','-append');

                A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
                AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
                option.mxitr=opt.maxItr;
                option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
                [s, out] = FPC_AS_mod(length(At(conf.y)),AO,conf.y,mu,[],option);
                fpcas{i,j,k}=out; fpcas{i,j,k}.alpha = conf.Psi(s);
                fpcas{i,j,k}.fVal(1)=0.5*sqrNorm(conf.Phi( fpcas{i,j,k}.alpha)-conf.y);
                fpcas{i,j,k}.fVal(2)=sqrNorm( fpcas{i,j,k}.alpha.*(fpcas{i,j,k}.alpha<0));
                fpcas{i,j,k}.fVal(3)=pNorm(conf.Psit(fpcas{i,j,k}.alpha),1);
                fpcas{i,j,k}.opt = opt; alphaHat=fpcas{i,j,k}.alpha;
                fpcas{i,j,k}.RMSE=sqrNorm(alphaHat-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
                fprintf('fpcas RMSE=%g\n',fpcas{i,j,k}.RMSE);
                save(filename,'fpcas','-append');
            end
        end
    end
end

if(any(runList==915))
    load(filename,'*015'); j=1;
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7];
    K = 1;

    npgTime   = 0;
    npgcTime  = 0;
    npgsTime  = 0;
    fistaTime = 0;
    spiralTime= 0;
    fpcasTime = 0;
    npgCost   = 0;
    npgcCost  = 0;
    npgsCost  = 0;
    fistaCost = 0;
    spiralCost= 0;
    fpcasCost = 0;
    npgRMSE   = 0;
    npgcRMSE  = 0;
    npgsRMSE  = 0;
    fistaRMSE = 0;
    spiralRMSE= 0;
    fpcasRMSE = 0;

    for k=1:K
        npgTime   =   npgTime+1/K*showResult(   npg015(:,:,k),2,'time');
        npgcTime  =  npgcTime+1/K*showResult(  npgc015(:,:,k),2,'time');
        npgsTime  =  npgsTime+1/K*showResult(  npgs015(:,:,k),2,'time');
        fistaTime = fistaTime+1/K*showResult( fista015(:,:,k),2,'time');
        spiralTime=spiralTime+1/K*showResult(spiral015(:,:,k),2,'time');
        fpcasTime = fpcasTime+1/K*showResult( fpcas015(:,:,k),2,'cpu');

        npgCost   =   npgCost+1/K*showResult(   npg015(:,:,k),2,'cost');
        npgcCost  =  npgcCost+1/K*showResult(  npgc015(:,:,k),2,'cost');
        npgsCost  =  npgsCost+1/K*showResult(  npgs015(:,:,k),2,'cost');
        fistaCost = fistaCost+1/K*showResult( fista015(:,:,k),2,'cost');
        spiralCost=spiralCost+1/K*showResult(spiral015(:,:,k),2,'cost');
        fpcasCost = fpcasCost+1/K*showResult( fpcas015(:,:,k),2,'f');

        npgRMSE   =   npgRMSE+1/K*showResult(   npg015(:,:,k),2,'RMSE');
        npgcRMSE  =  npgcRMSE+1/K*showResult(  npgc015(:,:,k),2,'RMSE');
        npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs015(:,:,k),2,'RMSE');
        fistaRMSE = fistaRMSE+1/K*showResult( fista015(:,:,k),2,'RMSE');
        spiralRMSE=spiralRMSE+1/K*showResult(spiral015(:,:,k),2,'reconerror');
        fpcasRMSE = fpcasRMSE+1/K*showResult( fpcas015(:,:,k),2,'RMSE');
    end

    % npgTime(:,1:2)=[];
    % npgcTime(:,1:2)=[];
    %  npgsTime(:,1:2)=[];
    % spiralTime(:,1:2)=[];
    % fpcasTime(:,1:2)=[];
    % npgCost(:,1:2)=[];
    % npgcCost(:,1:2)=[];
    %  npgsCost(:,1:2)=[];
    % spiralCost(:,1:2)=[];
    % fpcasCost(:,1:2)=[];
    % npgRMSE(:,1:2)=[];
    % npgcRMSE(:,1:2)=[];
    %  npgsRMSE(:,1:2)=[];
    % spiralRMSE(:,1:2)=[];
    % fpcasRMSE(:,1:2)=[];
    % fpcasRMSE_(:,1:2)=[];
    %  npgsRMSE_(:,1:2)=[];

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r); [r,c1(idx1)]
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r); [r,c2(idx2)]
    [r,c3]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r); [r,c3(idx3)]
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r); [r,c4(idx4)]
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r); [r,c5(idx5)]
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r); [r,c6(idx6)]

    figure;
    loglog(snr,   npgRMSE(:,3),'r-*'); hold on;
    loglog(snr,  npgcRMSE(:,3),'c-p');
    loglog(snr,  npgsRMSE(:,3),'k-s');
    loglog(snr,spiralRMSE(:,3),'k-^');
    loglog(snr, fpcasRMSE(:,3),'g-<');
    loglog(snr, fistaRMSE(:,3),'b-.');
    
    figure;
    loglog(snr,   npgTime(:,3),'r-*'); hold on;
    loglog(snr,  npgcTime(:,3),'c-p');
    loglog(snr,  npgsTime(:,3),'k-s');
    loglog(snr,spiralTime(:,3),'k-^');
    loglog(snr, fpcasTime(:,3),'g-<');
    loglog(snr, fistaTime(:,3),'b-.');
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
    opt.continuation=false;
    for k=1:1
        for i=7:length(m)
            opt.m=m(i); opt.snr=inf;
            opt=loadLinear(conf,opt);
            initSig = conf.Phit(conf.y)*0+1;
            fprintf('min=%d, max=%d\n',min(conf.y), max(conf.y));

            for j=1:5
                opt.u = u(i)*10^(j-3);
                fprintf('%s, i=%d, j=%d, k=%d\n','FISTA_ADMM_NNL1',i,j,k);

                % opt.alphaStep='FISTA_ADMM_NNL1';
                % npg{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'npg','-append');

                % opt.alphaStep='IST_ADMM_NNL1';
                % ist{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'ist','-append');




                opt.alphaStep='FISTA_L1';
                npgs{i,j,k}=lasso(conf.Phi,conf.Phit,...
                    conf.Psi,conf.Psit,conf.y,initSig,opt);
                save(filename,'npgs','-append');

                % opt.alphaStep='FISTA_L1';
                % opt.initStep='fixed'; opt.adaptiveStep=false;
                % fista{i,j,k}=lasso(conf.Phi,conf.Phit,...
                %     conf.Psi,conf.Psit,conf.y,initSig,opt);
                % save(filename,'fista','-append');
                % opt=rmfield(opt,'initStep'); opt=rmfield(opt,'adaptiveStep');

                % subtolerance=1e-5;
                % [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                %     SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                %     'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','poisson',...
                %     'initialization',initSig,'maxiter',opt.maxItr,...
                %     'miniter',0,'stopcriterion',3,...
                %     'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                %     'subtolerance',subtolerance,'monotone',1,...
                %     'saveobjective',1,'savereconerror',1,'savecputime',1,...
                %     'reconerrortype',3,'savedifalpha',1,...
                %     'savesolutionpath',0,'verbose',1000);
                % out.opt=opt; spiral{i,j,k}=out;
                % save(filename,'spiral','-append');
            end
        end
    end
end

if(any(runList==906))
    filename = [mfilename '_006.mat']; load(filename);

    K=5;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    fprintf('Poisson example\n');
    npgTime   = 0;
    npgsTime  = 0;
    fistaTime = 0;
    istTime   = 0;
    spiralTime= 0;
    npgCost   = 0;
    npgsCost  = 0;
    fistaCost = 0;
    istCost   = 0;
    spiralCost= 0;
    npgRMSE   = 0;
    npgsRMSE  = 0;
    fistaRMSE = 0;
    istRMSE   = 0;
    spiralRMSE= 0;

    for k=1:K
        npgTime   =   npgTime+1/K*showResult(   npg(:,:,k),2,'time');
        npgsTime  =  npgsTime+1/K*showResult(  npgs(:,:,k),2,'time');
        fistaTime = fistaTime+1/K*showResult( fista(:,:,k),2,'time');
        spiralTime=spiralTime+1/K*showResult(spiral(:,:,k),2,'time');
        istTime   =   istTime+1/K*showResult(   ist(:,:,k),2,'time');

        npgCost   =   npgCost+1/K*showResult(   npg(:,:,k),2,'cost');
        npgsCost  =  npgsCost+1/K*showResult(  npgs(:,:,k),2,'cost');
        fistaCost = fistaCost+1/K*showResult( fista(:,:,k),2,'cost');
        spiralCost=spiralCost+1/K*showResult(spiral(:,:,k),2,'cost');
        istCost   =   istCost+1/K*showResult(   ist(:,:,k),2,'cost');

        npgRMSE   =   npgRMSE+1/K*showResult(   npg(:,:,k),2,'RMSE');
        npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs(:,:,k),2,'RMSE');
        fistaRMSE = fistaRMSE+1/K*showResult( fista(:,:,k),2,'RMSE');
        spiralRMSE=spiralRMSE+1/K*showResult(spiral(:,:,k),2,'reconerror');
        istRMSE   =   istRMSE+1/K*showResult(   ist(:,:,k),2,'RMSE');

        % npgRMSE   =   npgRMSE+1/K*showResult(   npg(:,:,k),4,2);
        % npgsRMSE  =  npgsRMSE+1/K*showResult(  npgs(:,:,k),4,2);
        % fistaRMSE = fistaRMSE+1/K*showResult( fista(:,:,k),4,2);
        % spiralRMSE=spiralRMSE+1/K*showResult(spiral(:,:,k),4,2);
        % istRMSE   =   istRMSE+1/K*showResult(   ist(:,:,k),4,2);
    end

    [r,c1]=find(   npgRMSE==repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c3]=find( fistaRMSE==repmat(min( fistaRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c5]=find(   istRMSE==repmat(min(   istRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [c1(idx1) ,c2(idx2) ,c3(idx3) ,c4(idx4) ,c5(idx5)]
    idx1 = 3;    
    idx2 = 3;
    idx3 = 5;
    idx4 = 3;
    idx5 = 3;

    figure;
    semilogy(m,   npgRMSE(:,idx1),'r-*'); hold on;
    loglog(m,  npgsRMSE(:,idx2),'c-p');
    %loglog(m, fistaRMSE(:,idx3),'g-s');
    loglog(m,spiralRMSE(:,idx4),'k-^');
    %loglog(m,   istRMSE(:,idx5),'b-.');

    figure;
    plot(m,   npgTime(:,idx1),'r-*'); hold on;
    loglog(m,  npgsTime(:,idx2),'c-p');
    %loglog(m, fistaTime(:,idx3),'g-s');
    loglog(m,spiralTime(:,idx4),'k-^');
    %loglog(m,   istTime(:,idx5),'b-.');

    keyboard;

    forSave=[];
    forSave=[forSave,    npgTime(:,idx1)];
    forSave=[forSave,   npgsTime(:,idx2)];
    forSave=[forSave,  fistaTime(:,idx3)];
    forSave=[forSave, spiralTime(:,idx4)];
    forSave=[forSave,    istTime(:,idx5)];

    forSave=[forSave,    npgCost(:,idx1)];
    forSave=[forSave,   npgsCost(:,idx2)];
    forSave=[forSave,  fistaCost(:,idx3)];
    forSave=[forSave, spiralCost(:,idx4)];
    forSave=[forSave,    istCost(:,idx5)];

    forSave=[forSave,    npgRMSE(:,idx1)];
    forSave=[forSave,   npgsRMSE(:,idx2)];
    forSave=[forSave,  fistaRMSE(:,idx3)];
    forSave=[forSave, spiralRMSE(:,idx4)];
    forSave=[forSave,    istRMSE(:,idx5)];
    forSave=[forSave, m(:)];
    save('varyMeasurementPoisson.data','forSave','-ascii');

    mIdx=5; forSave=[]; t=0;
    t=t+1; temp=   npg{mIdx,idx1,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgs{mIdx,idx2,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= fista{mIdx,idx3,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=   ist{mIdx,idx5,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    save('stepSize.data','forSave','-ascii');

    forSave=[]; t=0;
    out=   ist{mIdx,idx5,1};
    t=t+1; forSave(1:length(out.cost      ),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE      ),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time      ),t)=out.time;
    out=   npg{mIdx,idx1,1};
    t=t+1; forSave(1:length(out.cost      ),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE      ),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time      ),t)=out.time;
    out=spiral{mIdx,idx4,1};
    t=t+1; forSave(1:length(out.cost      ),t)=out.cost;
    t=t+1; forSave(1:length(out.reconerror),t)=out.reconerror;
    t=t+1; forSave(1:length(out.time      ),t)=out.time;

    keyboard
    mincost=reshape(forSave(:,[1,4,7]),[],1); 
    mincost=min(mincost(mincost~=0));
    idx=(forSave(:,1)~=0); forSave(idx,1)=(forSave(idx,1)-mincost);
    idx=(forSave(:,4)~=0); forSave(idx,4)=(forSave(idx,4)-mincost);
    idx=(forSave(:,7)~=0); forSave(idx,7)=(forSave(idx,7)-mincost);


    save('cost_itr.data','forSave','-ascii');
    system(['mv varyMeasurementPoisson.data stepSize.data cost_itr.data ' paperDir]);
end

% The X-ray CT example, test and find the best u for each prjFull
if(any(runList==008))     % FPCAS
    filename = [mfilename '_008.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u=10.^[-6 -5 -5 -5 -5 -5];
    opt.maxItr=2e3; opt.thresh=1e-6;
    for i=1:length(prjFull)
        fprintf('%s, i=%d, j=%d\n','X-ray CT example',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2; opt.u = u(i);
        opt=conf.setup(opt);
        conf.y=conf.Phi(opt.trueAlpha); % equivalent to linear projection
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

        tic
        fbp{i,j}.img=conf.FBP(conf.y);
        fbp{i,j}.time=toc;
        fbp{i,j}.alpha=fbp{i,j}.img(opt.mask~=0);
        fbp{i,j}.RSE=sqrNorm(conf.y-conf.Phi(fbp{i,j}.alpha))/sqrNorm(conf.y);
        %fbp{i,j}.RMSE=1-(fbp{i,j}.alpha'*opt.trueAlpha/norm(fbp{i,j}.alpha)/norm(opt.trueAlpha))^2;
        fbp{i,j}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);
        save(filename,'fbp','-append');

        % opt.alphaStep='FISTA_ADMM_NNL1';
        % npg{i,j}=lasso(conf.Phi,conf.Phit,...
        %     conf.Psi,conf.Psit,conf.y,initSig,opt);
        % save(filename,'npg','-append');

        % A = @(xx) conf.Phi(conf.Psi(xx)); At = @(yy) conf.Psit(conf.Phit(yy));
        % AO=A_operator(A,At); mu=opt.u; option.x0=conf.Psit(initSig);
        % option.mxitr=opt.maxItr;
        % option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
        % [s, out] = FPC_AS_mod(length(At(conf.y)),AO,conf.y,mu,[],option);
        % fpcas{i,j}=out; fpcas{i,j}.alpha = conf.Psi(s);
        % fpcas{i,j}.fVal(1)=0.5*sqrNorm(conf.Phi(fpcas{i,j}.alpha)-conf.y);
        % fpcas{i,j}.fVal(2)=sqrNorm(fpcas{i,j}.alpha.*(fpcas{i,j}.alpha<0));
        % fpcas{i,j}.fVal(3)=pNorm(conf.Psit(fpcas{i,j}.alpha),1);
        % fpcas{i,j}.opt = opt;
        % fpcas{i,j}.RMSE=sqrNorm(fpcas{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        % fprintf('fpcas RMSE=%g\n',fpcas{i}.RMSE);
        % save(filename,'fpcas','-append');

        % out=[]; subtolerance=1e-5;
        % [out.alpha, out.p, out.cost, out.reconerror, out.time, out.difAlpha] = ...
        %     SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
        %     'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
        %     'initialization',initSig,'maxiter',opt.maxItr,...
        %     'miniter',0,'stopcriterion',3,...
        %     'tolerance',opt.thresh,'truth',opt.trueAlpha,...
        %     'subtolerance',subtolerance,'monotone',1,...
        %     'saveobjective',1,'savereconerror',1,'savecputime',1,...
        %     'reconerrortype',3,'savedifalpha',1,...
        %     'savesolutionpath',0,'verbose',10);
        % out.opt=opt; spiral{i,j}=out;
        % spiral{i,j}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral{i,j}.alpha)-conf.y);
        % spiral{i,j}.fVal(2)=sqrNorm(spiral{i,j}.alpha.*(spiral{i,j}.alpha<0));
        % spiral{i,j}.fVal(3)=sum(abs(conf.Psit(spiral{i,j}.alpha)));
        % save(filename,'spiral','-append');
    end
end

if(any(runList==908))
    filename = [mfilename '_008.mat']; load(filename);
    i=3; t=0;
    prjFull = [60, 80, 100, 120, 180, 360];
    perform=[]; forSave=[];

    out=   npg{i};
    fprintf('NPGwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'NPGwLin.png','png');
    temp=showResult(   npg,2,'RMSE');
    perform=[perform, temp];
    %perform=[perform, [temp(1,6); temp(2:end,5)]];
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    out=spiral{i};
    fprintf('SPIRALwLin: %g\n',out.reconerror(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'SPIRALwLin.png','png');
    temp=showResult(spiral,2,'reconerror');
    perform=[perform, temp(:)];
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.reconerror),t)=out.reconerror;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    out=fbp{i};
    fprintf('FBPwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,spiral{i}.opt.mask);
    imwrite(img/max(img(:)),'FBPwLin.png','png');
    temp=showResult(fbp,2,'RMSE');
    perform=[perform, temp(:)];

    out=fpcas{i};
    fprintf('FPCASwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'FPCASwLin.png','png');
    temp=showResult(fpcas,2,'RMSE');
    perform=[perform, temp(:)];

    temp=showResult(   npg,2,'time'); perform=[perform, temp(:)];
    temp=showResult(spiral,2,'time'); perform=[perform, temp(:)];
    temp=showResult(   fbp,2,'time'); perform=[perform, temp(:)];
    temp=showResult( fpcas,2, 'cpu'); perform=[perform, temp(:)];

    figure;
    for i=1:4
        semilogy(prjFull/2,perform(:,i)); hold on;
    end

    keyboard

    perform=[perform,prjFull(:)];
    save('varyPrj.data','perform','-ascii');
    save('xray_itr.data','forSave','-ascii');
    !cp varyPrj.data xray_itr.data *wLin.png ~/research/myPaper/asilomar2014/
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
        conf.y = conf.y/max(conf.y(:)); % important to normalize
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


if(any(runList==909))
    filename = [mfilename '_009.mat']; load(filename);

    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    fprintf('Poisson Log link example with glass beads\n');

    k=1;
    npgTime   = showResult(   npg(:,:,k),2,'time');
    npgsTime  = showResult(  npgs(:,:,k),2,'time');
    % fbpTime   = showResult(   fbp(:,:,k),2,'time');

    npgCost   = showResult(   npg(:,:,k),2,'cost');
    npgsCost  = showResult(  npgs(:,:,k),2,'cost');
    % fbpCost   = showResult(   fbp(:,:,k),2,'cost');

    npgRMSE   = showResult(   npg(:,:,k),2,'RMSE');
    npgsRMSE  = showResult(  npgs(:,:,k),2,'RMSE');
    fbpRMSE   = showResult(   fbp(:,:,k),2,'RMSE');

    % [r,c1]=find(   npgRMSE==repmat(min(   npgRMSE,[],2),1,3)); [r,idx1]=sort(r);
    % [r,c2]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,3)); [r,idx2]=sort(r);
    % [r,c3]=find(   fbpRMSE==repmat(min(   fbpRMSE,[],2),1,3)); [r,idx3]=sort(r);
    % [c1(idx1) ,c2(idx2) ,c3(idx3)]
    idx1 = 2;    
    idx2 = 2;
    idx3 = 2;

    figure;
    semilogy(prjFull,   npgRMSE(:,idx1),'r-*'); hold on;
    loglog(prjFull,  npgsRMSE(:,idx2),'c-p');
    loglog(prjFull, fbpRMSE(:,idx3),'g-s');

    figure;
    plot(prjFull,   npgTime(:,idx1),'r-*'); hold on;
    loglog(prjFull,  npgsTime(:,idx2),'c-p');
    %loglog(prjFull,   fbpTime(:,idx3),'g-s');

    forSave=[];
    forSave=[forSave,    npgRMSE(:,idx1)];
    forSave=[forSave,   npgsRMSE(:,idx2)];
    forSave=[forSave,    fbpRMSE(:,idx3)];

    forSave=[forSave,    npgTime(:,idx1)];
    forSave=[forSave,   npgsTime(:,idx2)];
    % forSave=[forSave,    fbpTime(:,idx3)];

    % forSave=[forSave,    npgCost(:,idx1)];
    % forSave=[forSave,   npgsCost(:,idx2)];
    % forSave=[forSave,    fbpCost(:,idx3)];

    forSave=[forSave, prjFull(:)];
    save('varyPrjGlassBead.data','forSave','-ascii');

    idx=4;
    img=showImgMask(npg {idx,idx1}.alpha,npg{idx,idx1}.opt.mask);
    imwrite(img/max(img(:)),'NPGgb.png','png');
    img=showImgMask(npgs{idx,idx2}.alpha,npg{idx,idx1}.opt.mask);
    imwrite(img/max(img(:)),'NPGSgb.png','png');
    img=showImgMask(fbp {idx,idx3}.alpha,npg{idx,idx1}.opt.mask);
    imwrite(img/max(img(:)),'FBPgb.png','png');

    disp([npg{idx,idx1}.RMSE(end), npgs{idx,idx2}.RMSE(end), fbp{idx,idx3}.RMSE]);

    system(['mv varyPrjGlassBead.data NPGgb.png NPGSgb.png FBPgb.png ' paperDir]);
    keyboard
end

end

function e=gEle(x,i); e=x(i); end

