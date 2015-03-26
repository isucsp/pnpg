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


end

