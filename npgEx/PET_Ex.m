function PET_Ex(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%
%
%                  PET example, with background noise b
%      Vary the total counts of the measurements, with continuation

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        % PET example
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt'); filename = [mfilename '.mat'];
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1; opt.noiseType='poisson';

        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        K=3;

        as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
        a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];
        atv= [-0.5,-0.5,  0,   0, 0.5,   1];

        opt.mask=[];
        for k=1:K
            for i=1:length(count)
                [y,Phi,Phit,Psi,Psit,fbpfunc,opt]=loadPET(count(i),opt);

                fbp{i,1,k}.alpha=maskFunc(fbpfunc(y),opt.mask~=0);
                fbp{i,1,k}.RMSE=sqrNorm(fbp{i,1,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);

                fprintf('fbp RMSE=%f\n',sqrNorm(fbp{i,1,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha));
                fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
                u_max=1;

                initSig=max(fbp{i,1,k}.alpha,0);

                opt.fullcont=false;
                j=1;
                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);
                opt.contShrnk=0.1; opt.contGamma=15;

                opt.u = 10^atv(i)*u_max; opt.proximal='tviso';
                if(k==1 && i==5)
                    Opt=opt; opt.adaptiveStep=false;
                    % npgTV_noAdpStp=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt.thresh=1e-10;
                    %npgTV_noAdpStpLong=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt=Opt; opt.cumuTol=1; opt.incCumuTol=false;
                    %npgTV_n1=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt.incCumuTol=true;
                    test=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    keyboard
                    opt=Opt;
                    npgTV_n4=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt.thresh=1e-10;
                    spiralTV_Long=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);



                save(filename);
                    keyboard
                else
                    continue
                end

%               npgTV {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgTVc{i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               spiralTV{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);

                save(filename);

                opt.u = 10^a(i)*u_max; opt.proximal='wvltADMM';
%               npg   {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgc  {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.proximal='wvltLagrangian';
%               spiral{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);

%               opt.u = 10^as(i)*u_max; opt.proximal='wvltADMM';
%               npgs  {i,j,k}=Wrapper.NPGs   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgsc {i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);

                save(filename);

%               opt.fullcont=true;
%               % for isotv
%               u_max=1;
%               aa =(3:-0.5:-6);
%               opt.u=(10.^aa)*u_max; opt.proximal='tviso';
%               if(i<5) continue; end
%               npgTVFull{i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt);
%               for j=1:length(aa); if(aa(j)>-2)
%                   opt.u=10^aa(j)*u_max; opt.proximal='tviso';
%                   if(j==1)
%                       spiralTVFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
%                   else
%                       spiralTVFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,spiralTVFull{i,j-1,k}.alpha,opt);
%                   end
%               end; end

%               % for wavelet l1 norm
%               u_max=1;
%               aa = (3:-0.5:-6);
%               opt.u=(10.^aa)*u_max; opt.proximal='wvltADMM';
%               npgFull {i,k}=Wrapper.NPG (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgsFull{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt);
%               for j=1:length(aa); if(aa(j)>-2)
%                   opt.u=10^aa(j)*u_max; opt.proximal='wvltLagrangian';
%                   if(j==1)
%                       spiralFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
%                   else
%                       spiralFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,spiralFull{i,j-1,k}.alpha,opt);
%                   end
%               end; end

%               % following are methods for weighted versions
%               ty=max(sqrt(y),1);
%               wPhi=@(xx) Phi(xx)./ty;
%               wPhit=@(xx) Phit(xx./ty);
%               wy=(y-opt.bb(:))./ty;
%               wu_max=pNorm(Psit(wPhit(wy)),inf);
%               opt.noiseType='gaussian';

%               opt.fullcont=true;
%               opt.u=(10.^aa)*wu_max; opt.maxItr=1e4; opt.thresh=1e-12;
%               wnpgFull {i,k}=Wrapper.NPG(wPhi,wPhit,Psi,Psit,wy,initSig,opt); out=wnpgFull{i,k};
%               fprintf('k=%d, good a = 1e%g\n',k,max((aa(out.contRMSE==min(out.contRMSE)))));
%               opt.fullcont=false;

%               opt.u = 10^a(i)*u_max;
%               fprintf('%s, i=%d, j=%d, k=%d\n','PET Example_003',i,1,k);
%               wnpg{i,k}=Wrapper.NPG         (wPhi,wPhit,Psi,Psit,wy,initSig,opt);
%               wspiral{i,k}=Wrapper.SPIRAL (wPhi,wPhit,Psi,Psit,wy,initSig,opt);
%               % wnpgc  {i,k}=Wrapper.NPGc   (wPhi,wPhit,Psi,Psit,wy,initSig,opt);
            end
        end
    case 'bound'
        clear('opt');
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1; opt.noiseType='poisson';
        opt.minItr=30;
        opt.mask  =[];

        K=1;
        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        for k=1:K
            for i=1:length(count)
                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,1,k);
                [y,Phi,Phit,Psi,Psit,fbpfunc,opt]=loadPET(count(i),opt);

                [x0s,g]=Utils.poissonModelConstEst(Phi,Phit,y,opt.bb);
                u_maxANI(i,k)=TV.upperBoundU(maskFunc(g,opt.mask));
                u_maxISO(i,k)=sqrt(2)*u_maxANI(i,k);
                initSig=ones(size(opt.trueAlpha))*x0s;

                rmse=0; opt.u=u_maxANI(i,k); opt.proximal='tvl1';
                while(rmse==0)
                    opt.u = 0.9*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                opt.u=opt.u/0.9; rmse=0;
                while(rmse==0)
                    opt.u = 0.99*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                u_trueANI(i,k)=opt.u/0.99;

                rmse=0; opt.u=u_maxISO(i,k); opt.proximal='tviso';
                while(rmse==0)
                    opt.u = 0.9*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                opt.u=opt.u/0.9; rmse=0;
                while(rmse==0)
                    opt.u = 0.99*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                u_trueISO(i,k)=opt.u/0.99;
            end
        end

        figure;
        loglog(count,u_maxANI,'b^-'); hold on;
        loglog(count,u_trueANI,'bs--');
        loglog(count,u_maxISO,'r*-');
        loglog(count,u_trueISO,'ro--');
        h=legend('U_0','empirical anisotropic U','sqrt(2)U_0','empirical isotropic U');
        set(h,'interpreter','latex');

        forSave=[count(:) u_maxANI, u_trueANI, u_maxISO, u_trueISO];
        save('bound4U.data','forSave','-ascii');

        save('PET_Ex_bound.mat');

    case lower('plotTV')
        filename = [mfilename '.mat']; load(filename);
        fprintf('PET Poisson TV example\n');

        count = [1e4 1e5 1e6 1e7 1e8 1e9];

        K = 1:3;

        npgTVcTime= mean(Cell.getField(npgTVc(:,1,K),'time'),3);
        npgTVcCost= mean(Cell.getField(npgTVc(:,1,K),'cost'),3);
        npgTVcRMSE= mean(Cell.getField(npgTVc(:,1,K),'RMSE'),3);
        npgTVTime = mean(Cell.getField(npgTV (:,1,K),'time'),3);
        npgTVCost = mean(Cell.getField(npgTV (:,1,K),'cost'),3);
        npgTVRMSE = mean(Cell.getField(npgTV (:,1,K),'RMSE'),3);
        spiralTVTime = mean(Cell.getField(spiralTV (:,1,K),'time'),3);
        spiralTVCost = mean(Cell.getField(spiralTV (:,1,K),'cost'),3);
        spiralTVRMSE = mean(Cell.getField(spiralTV (:,1,K),'RMSE'),3);

        npgsTime  = mean(Cell.getField(  npgs(:,1,K),'time'),3);
        npgsCost  = mean(Cell.getField(  npgs(:,1,K),'cost'),3);
        npgsRMSE  = mean(Cell.getField(  npgs(:,1,K),'RMSE'),3);

        npgscTime = mean(Cell.getField( npgsc(:,1,K),'time'),3);
        npgscCost = mean(Cell.getField( npgsc(:,1,K),'cost'),3);
        npgscRMSE = mean(Cell.getField( npgsc(:,1,K),'RMSE'),3);

        fbpRMSE   = mean(Cell.getField(   fbp(:,1,K),'RMSE'),3);

        figure;
        loglog(count,npgTVRMSE,'r-*'); hold on;
        loglog(count,   fbpRMSE,'b-o');
        loglog(count,spiralTVRMSE,'k-^');
        loglog(count,  npgTVcRMSE,'k*-.');
        loglog(count,  npgsRMSE,'c>-');
        loglog(count, npgscRMSE,'gs-');
        legend('npgTV','fbp','spiralTV','npgTVc','npgs','npgsc');

        figure;
        loglog(count,   npgTVTime,'r-*'); hold on;
        loglog(count,spiralTVTime,'k-^');
        loglog(count,  npgTVcTime,'k*-.');
        loglog(count,  npgsTime,'c>-');
        loglog(count, npgscTime,'gs-');
        legend('npgTV','spiralTV','npgTVc','npgs','npgsc');

        forSave=[npgTVTime, npgTVcTime, npgsTime, npgscTime, spiralTVTime,...
            npgTVCost, npgTVcCost, npgsCost, npgscCost, spiralTVCost,...
            npgTVRMSE, npgTVcRMSE, npgsRMSE, npgscRMSE, spiralTVRMSE,...
            fbpRMSE, count(:)];
        save('varyCntPETTV.data','forSave','-ascii');

        forSave=[]; t=0; mIdx=5; k=1;
        out=   npgTV_n4;
        t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=spiralTV_Long;
        t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=npgTV_noAdpStpLong;
        t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=   npgTV_n1;
        t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        save('stepSize.data','forSave','-ascii');

        figure; semilogy(forSave(:,1),'r'); hold on;
        semilogy(forSave(:,3),'g');
        semilogy(forSave(:,5),'b');
        title('step size versus number of iterations');
        legend('npgTV','spiralTV','npgTV noAdaptive Step');

        forSave=[]; t=0; mIdx=5; k=1;
        out=   npgTV_n4;
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=  npgTV_n1;
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=  npgTV_noAdpStpLong;
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        t=t+1; forSave(1:length(out.difAlpha),t)=out.difAlpha;
        out=  spiralTV_Long;
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        t=t+1; forSave(1:length(out.difAlpha),t)=out.difAlpha;

        save('cost_itrPETTV.data','forSave','-ascii');
        mincost=reshape(forSave(:,[1,4,7,11]),[],1); 
        mincost=min(mincost(mincost~=0));

        figure;
        semilogy(forSave(:,3),forSave(:,1)-mincost,'r'); hold on;
        semilogy(forSave(:,6),forSave(:,4)-mincost,'g');
        semilogy(forSave(:,9),forSave(:,7)-mincost,'b');
        semilogy(forSave(:,13),forSave(:,11)-mincost,'k');
        legend('npgTV n4','npgTV n1','npgTV nInf','spiralTV');
        hold on;
        idx=min(find(forSave(:,10)<1e-6));
        plot(forSave(idx,9),forSave(idx,7)-mincost,'bo');
        xxx=idx;
        idx=min(find(forSave(10:end,14)<1e-6))+10;
        plot(forSave(idx,13),forSave(idx,11)-mincost,'k*');
        xxx=[xxx;idx];  xxx=xxx(:)';
        save('cost_itrPETTVidx.data','xxx','-ascii');

        figure; semilogy(forSave(:,3),forSave(:,2),'r'); hold on;
        semilogy(forSave(:,6),forSave(:,5),'g');
        semilogy(forSave(:,9),forSave(:,8),'b');
        semilogy(forSave(:,13),forSave(:,12),'k');
        legend('npgTV n4','npgTV n1','npgTV nInf','spiralTV');

        keyboard

        nn=128;
        xtrue = read_zubal_emis('nx', nn, 'ny', nn);
        % attenuation map
        mumap = read_zubal_attn('nx', nn, 'ny', nn);
        imwrite(xtrue/max(xtrue(:)),'pet.png');
        imwrite(mumap/max(mumap(:)),'mumap.png');

        idx=5;
        fprintf('   NPGTV: %g%%\n',   npgTV{idx}.RMSE(end)*100);
        fprintf('SPIRALTV: %g%%\n',spiralTV{idx}.RMSE(end)*100);
        fprintf('     FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueAlpha)*100);
        fprintf('    NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
        img=npg{idx}.alpha; mask=npg{idx}.opt.mask;
        img=showImgMask(   npgTV{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPGTV_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPGTV_pet.png')
        img=showImgMask(spiralTV{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRALTV_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRALTV_pet.png')
        img=showImgMask(     fbp{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,     'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),     'FBP_pet.png')
        img=showImgMask(    npgs{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,    'NPGs_pet.eps','psc2'); imwrite(img/max(xtrue(:)),    'NPGs_pet.png')

        idx=4;
        fprintf('   NPGTV: %g%%\n',   npgTV{idx}.RMSE(end)*100);
        fprintf('SPIRALTV: %g%%\n',spiralTV{idx}.RMSE(end)*100);
        fprintf('     FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueAlpha)*100);
        fprintf('    NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
        img=npg{idx}.alpha; mask=npg{idx}.opt.mask;
        img=showImgMask(   npgTV{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPGTV_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPGTV_pet2.png')
        img=showImgMask(spiralTV{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRALTV_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRALTV_pet2.png')
        img=showImgMask(     fbp{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,     'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),     'FBP_pet2.png')
        img=showImgMask(    npgs{idx}.alpha,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,    'NPGs_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),    'NPGs_pet2.png')

        paperDir='~/research/myPaper/asilomar2015/';
        decide=input(sprintf('start to copy to %s [y/N]?',paperDir));
        if strcmpi(decide,'y')
            system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
        end
        system('rm *_pet.eps *_pet2.eps *_pet2.png');
        close all;

    case 'plot'
        filename = [mfilename '.mat'];
        load(filename);
        fprintf('PET Poisson example\n');

        count = [1e4 1e5 1e6 1e7 1e8 1e9];

        K = 1:1;

        npgTime      = mean(Cell.getField(npg     (:,1,K),'time'),3);
        npgCost      = mean(Cell.getField(npg     (:,1,K),'cost'),3);
        npgRMSE      = mean(Cell.getField(npg     (:,1,K),'RMSE'),3);
        npgcTime     = mean(Cell.getField(npgc    (:,1,K),'time'),3);
        npgcCost     = mean(Cell.getField(npgc    (:,1,K),'cost'),3);
        npgcRMSE     = mean(Cell.getField(npgc    (:,1,K),'RMSE'),3);
        npgsTime     = mean(Cell.getField(npgs    (:,1,K),'time'),3);
        npgsCost     = mean(Cell.getField(npgs    (:,1,K),'cost'),3);
        npgsRMSE     = mean(Cell.getField(npgs    (:,1,K),'RMSE'),3);
        npgscTime    = mean(Cell.getField(npgsc   (:,1,K),'time'),3);
        npgscCost    = mean(Cell.getField(npgsc   (:,1,K),'cost'),3);
        npgscRMSE    = mean(Cell.getField(npgsc   (:,1,K),'RMSE'),3);
        spiralTime   = mean(Cell.getField(spiral  (:,1,K),'time'),3);
        spiralCost   = mean(Cell.getField(spiral  (:,1,K),'cost'),3);
        spiralRMSE   = mean(Cell.getField(spiral  (:,1,K),'RMSE'),3);

        fbpRMSE      = mean(Cell.getField(fbp     (:,1,K),'RMSE'),3);

        figure;
        loglog(count,     npgRMSE,'r-*'); hold on;
        loglog(count,     fbpRMSE,'b-o');
        loglog(count,  spiralRMSE,'k-^');
        loglog(count,    npgcRMSE,'k*-.');
        loglog(count,    npgsRMSE,'c>-');
        loglog(count,   npgscRMSE,'gs-');
        legend('npg','fbp','spiral','npgc','npgs','npgsc');

        figure;
        loglog(count,     npgTime,'r-*'); hold on;
        loglog(count,  spiralTime,'k-^');
        loglog(count ,   npgcTime,'k*-.');
        loglog(count,    npgsTime,'c>-');
        loglog(count,   npgscTime,'gs-');
        legend('npg','spiral','npgc','npgs','npgsc');

        forSave=[npgTime, npgcTime, npgsTime, npgscTime, spiralTime,...
            npgCost, npgcCost, npgsCost, npgscCost, spiralCost,...
            npgRMSE, npgcRMSE, npgsRMSE, npgscRMSE, spiralRMSE,...
            fbpRMSE, count(:)...
            ];
        save('varyCntPET.data','forSave','-ascii');

        forSave=[]; t=0; mIdx=5; k=1;
        out=  npgc{mIdx,1,k};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=   npg{mIdx,1,k};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=spiral{mIdx,1,k};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=npgs{mIdx,1,k};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=npgsc_s;
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;

        save('cost_itrPET.data','forSave','-ascii');
        mincost=reshape(forSave(:,[1,4,7]),[],1); 
        mincost=min(mincost(mincost~=0));

        figure;
        semilogy(forSave(:,3),forSave(:,1)-mincost,'r'); hold on;
        semilogy(forSave(:,6),forSave(:,4)-mincost,'g');
        semilogy(forSave(:,9),forSave(:,7)-mincost,'b');
        if(mIdx==5 && k==1)
            semilogy(forSave(:,15),forSave(:,13)-min(max(forSave(:,13),0)),'c:');
            legend('npgc','npg','spiral','npgsc');
        else
            legend('npgc','npg','spiral');
        end
        figure; semilogy(forSave(:,3),forSave(:,2),'r'); hold on;
        semilogy(forSave(:,6),forSave(:,5),'g'); semilogy(forSave(:,9),forSave(:,8),'b');
        legend('npgc','npg','spiral');

        keyboard

        nn=128;
        xtrue = read_zubal_emis('nx', nn, 'ny', nn);
        % attenuation map
        mumap = read_zubal_attn('nx', nn, 'ny', nn);
        imwrite(xtrue/max(xtrue(:)),'pet.png');
        imwrite(mumap/max(mumap(:)),'mumap.png');

        idx=5;
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
        paperDir='~/research/myPaper/asilomar2014/'
        system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
        system('rm *_pet.eps *_pet2.eps *_pet2.png');
        close all;
    case 'fullplot'
        filename = [mfilename '.mat'];
        load(filename);

        k=1;
        aa =(3:-0.5:-6);
        for i=1:length(count)
            npgContRMSE  {i,k} = [  npgFull{i,k}.contRMSE(:);  npgFull{i,k}.RMSE(end)]; out=npgContRMSE{i,k};
            fprintf('i=%d, good a = 1e%g NPG\n',i,max((aa(out==min(out)))));
            npgsContRMSE {i,k} = [ npgsFull{i,k}.contRMSE(:); npgsFull{i,k}.RMSE(end)]; out=npgsContRMSE{i,k};
            fprintf('i=%d, good a = 1e%g NPGs\n',i,max((aa(out==min(out)))));
            npgTVContRMSE{i,k} = [npgTVFull{i,k}.contRMSE(:);npgTVFull{i,k}.RMSE(end)]; out=npgTVContRMSE{i,k};
            fprintf('i=%d, good a = 1e%g NPG_TV\n',i,max((aa(out==min(out)))));
        end

        for i=1:length(count)
            figure;
            semilogy(aa(1:length(npgContRMSE{i})),npgContRMSE{i},'r-*'); hold on;
            semilogy(aa(1:length(npgsContRMSE{i})),npgsContRMSE{i},'g-o');
            semilogy(aa(1:length(npgTVContRMSE{i})),npgTVContRMSE{i},'b-s');
            title(num2str(i));
            legend('NPG','NPGs','NPG-TV');
            aaa(i)=min(npgContRMSE{i});
            bbb(i)=min(npgsContRMSE{i});
            ccc(i)=min(npgTVContRMSE{i});
        end
        figure; semilogy(aaa,'r-*'); hold on;
        semilogy(bbb,'g-o');
        semilogy(ccc,'b-s');
        title('rmse vs count');
        legend('NPG','NPGs','NPG-TV');


        keyboard

        K=1;
        m = [1e4 1e5 1e6 1e7 1e8 1e9];
        fprintf('Poisson example\n');

        npgTime    = mean(Cell.getField(    npg(:,:,1:K),'time'),3);
        npgcTime   = mean(Cell.getField(   npgc(:,:,1:K),'time'),3);
        npgsTime   = mean(Cell.getField(   npgs(:,:,1:K),'time'),3);
        npgscTime  = mean(Cell.getField(  npgsc(:,:,1:K),'time'),3);
        spiralTime = mean(Cell.getField( spiral(:,:,1:K),'time'),3);

        npgCost    = mean(Cell.getField(    npg(:,:,1:K),'cost'),3);
        npgcCost   = mean(Cell.getField(   npgc(:,:,1:K),'cost'),3);
        npgsCost   = mean(Cell.getField(   npgs(:,:,1:K),'cost'),3);
        npgscCost  = mean(Cell.getField(  npgsc(:,:,1:K),'cost'),3);
        spiralCost = mean(Cell.getField( spiral(:,:,1:K),'cost'),3);

        npgRMSE    = mean(Cell.getField(    npg(:,:,1:K),'RMSE'),3);
        npgcRMSE   = mean(Cell.getField(   npgc(:,:,1:K),'RMSE'),3);
        npgsRMSE   = mean(Cell.getField(   npgs(:,:,1:K),'RMSE'),3);
        npgscRMSE  = mean(Cell.getField(  npgsc(:,:,1:K),'RMSE'),3);
        spiralRMSE = mean(Cell.getField( spiral(:,:,1:K),'RMSE'),3);

        aIdx=4;
        figure;
        loglog(m,    npgRMSE(:,aIdx),'r-*'); hold on;
        loglog(m,   npgsRMSE(:,aIdx),'c-p');
        loglog(m, spiralRMSE(:,aIdx),'k-^');
        loglog(m,   npgcRMSE(:,aIdx),'k*-.');
        loglog(m,  npgscRMSE(:,aIdx),'bs-.');
        legend('npg','npgs','spiral','npgc','npgsc');

        figure;
        loglog(m,    npgTime(:,aIdx),'r-*' ); hold on;
        loglog(m,   npgsTime(:,aIdx),'c-p' );
        loglog(m, spiralTime(:,aIdx),'k-^' );
        loglog(m,   npgcTime(:,aIdx),'k*-.');
        loglog(m,  npgscTime(:,aIdx),'bs-.');
        legend('npg','npgs','spiral','npgc','npgsc');

        forSave=[npgTime(:,aIdx), npgsTime(:,aIdx), npgcTime(:,aIdx), npgscTime(:,aIdx), spiralTime(:,aIdx),...
            npgCost(:,aIdx), npgsCost(:,aIdx), npgcCost(:,aIdx), npgscCost(:,aIdx), spiralCost(:,aIdx),...
            npgRMSE(:,aIdx), npgsRMSE(:,aIdx), npgcRMSE(:,aIdx), npgscRMSE(:,aIdx), spiralRMSE(:,aIdx),...
            m(:)];
        save('varyMeasurementPoisson.data','forSave','-ascii');


end

end
