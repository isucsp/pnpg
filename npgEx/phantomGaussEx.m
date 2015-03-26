function phantomGaussEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration
%
%                    Skyline log link Poisson Example
%    vary the number of measurements and inital intensity constatn I_0

% Linear model with Gaussian noise for Phantom example
% vary the number of measurements, and noise variance

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt'); filename = [mfilename '.mat'];
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

    case 'plot'
        load([mfilename '.mat']);

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

end

