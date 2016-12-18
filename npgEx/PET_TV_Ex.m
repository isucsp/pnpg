function PET_TV_Ex(op)
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
        clear -regexp '(?i)opt'
        filename = [mfilename '.mat'];
        OPT.mask=[]; OPT.outLevel=1;
        OPT.maxItr=1e4; OPT.thresh=1e-6; OPT.debugLevel=2; OPT.noiseType='poisson';
        %OPT.maxItr=10;

        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        K=1;

        as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
        a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];
        atv= [-0.5,-0.5,  0,   0, 0.5,   1];

        for k=1:K
            for i=length(count):-1:1
                j=1;
                [y,Phi,Phit,~,~,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
                NLL=@(x) Utils.poissonModel(x,Phi,Phit,y,OPT.bb);
                proximal=sparseProximal('iso',@(x)max(0,x));

                fbp{i,1,k}.x=fbpfunc(y);
                fbp{i,1,k}.RMSE=sqrNorm(fbp{i,1,k}.x-OPT.trueX)/sqrNorm(OPT.trueX);

                fprintf('fbp RMSE=%f\n',sqrNorm(fbp{i,1,k}.x-OPT.trueX)/sqrNorm(OPT.trueX));
                fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
                u_max=1;
                OPT.u = 10^atv(i)*u_max; OPT.proximal='tviso';
                OPT.stepShrnk=0.5; OPT.stepIncre=0.5;

                initSig=max(fbp{i,1,k}.x,0);

                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);
                if(k==1 && any(i==[4 6]))
                    opt=OPT;
                    pnpg_   {i,j,k}=pnpg(NLL,proximal,initSig,opt);
                end

                    continue;
                if(k==1 && any(i==[4 6]))
                    opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-5;
                    tfocs_200_m5 {i,j,k}=Wrapper.tfocs    (Phi,Phit,[],[],y,initSig,opt);
                    opt=OPT; opt.innerThresh=1e-5;
                    spiral_m5 {i,j,k}=Wrapper.SPIRAL  (Phi,Phit,[],[],y,initSig,opt);
                    mysave;
                end
                if(k==1 && any(i==[4 6]))
                    opt=OPT;
                    pnpg_   {i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.gamma=5; opt.b=0;
                    pnpgG5A0{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.gamma=5; opt.b=1/4;
                    pnpgG5Aq{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.gamma=15; opt.b=0;
                    pnpgGfA0{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.gamma=15; opt.b=1/4;
                    pnpgGfAq{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.b=0;
                    pnpgA0  {i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.cumuTol=0;
                    pnpg_n0 {i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
                    pnpg_n0m0{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    opt=OPT; opt.adaptiveStep=false;
                    pnpg_nInf{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                    mysave;
                end
                continue

                if(k==1) continue; end
                opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-5;
                tfocs_200_m5 {i,j,k}=Wrapper.tfocs    (Phi,Phit,[],[],y,initSig,opt);

                mysave;
                continue


                opt=OPT;
                spiralTV{i,j,k}=Wrapper.SPIRAL (Phi,Phit,[],[],y,initSig,opt);

                opt=OPT; opt.adaptiveStep=false; opt.thresh=1e-10;
                pnpgTV_noAdpStpLong{i,j,k}=pnpg(NLL,proximal,initSig,opt);
                opt=OPT; opt.thresh=1e-10;
                spiralTV_Long=Wrapper.SPIRAL (Phi,Phit,[],[],y,initSig,opt);
                opt=OPT; opt.innerThresh=1e-6;
                spiral_m6 {i,j,k}=Wrapper.SPIRAL  (Phi,Phit,[],[],y,initSig,opt);
                opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-4;
                tfocs_200_m4 {i,j,k}=Wrapper.tfocs    (Phi,Phit,[],[],y,initSig,opt);
                opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-6;
                tfocs_200_m6 {i,j,k}=Wrapper.tfocs    (Phi,Phit,[],[],y,initSig,opt);


%               npgTV {i,j,k}=Wrapper.NPG    (Phi,Phit,[],[],y,initSig,opt);
%               npgTVc{i,j,k}=Wrapper.NPGc   (Phi,Phit,[],[],y,initSig,opt);

                mysave;

%               opt.fullcont=true;
%               % for isotv
%               u_max=1;
%               aa =(3:-0.5:-6);
%               opt.u=(10.^aa)*u_max; opt.proximal='tviso';
%               if(i<5) continue; end
%               npgTVFull{i,k}=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
%               for j=1:length(aa); if(aa(j)>-2)
%                   opt.u=10^aa(j)*u_max; opt.proximal='tviso';
%                   if(j==1)
%                       spiralTVFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,[],[],y,initSig,opt);
%                   else
%                       spiralTVFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,[],[],y,spiralTVFull{i,j-1,k}.x,opt);
%                   end
%               end; end

%               % following are methods for weighted versions
%               ty=max(sqrt(y),1);
%               wPhi=@(xx) Phi(xx)./ty;
%               wPhit=@(xx) Phit(xx./ty);
%               wy=(y-opt.bb(:))./ty;
%               wu_max=pNorm([](wPhit(wy)),inf);
%               opt.noiseType='gaussian';

%               opt.fullcont=true;
%               opt.u=(10.^aa)*wu_max; opt.maxItr=1e4; opt.thresh=1e-12;
%               wnpgFull {i,k}=Wrapper.NPG(wPhi,wPhit,[],[],wy,initSig,opt); out=wnpgFull{i,k};
%               fprintf('k=%d, good a = 1e%g\n',k,max((aa(out.contRMSE==min(out.contRMSE)))));
%               opt.fullcont=false;

%               opt.u = 10^a(i)*u_max;
%               fprintf('%s, i=%d, j=%d, k=%d\n','PET Example_003',i,1,k);
%               wnpg{i,k}=Wrapper.NPG         (wPhi,wPhit,[],[],wy,initSig,opt);
%               wspiral{i,k}=Wrapper.SPIRAL (wPhi,wPhit,[],[],wy,initSig,opt);
%               % wnpgc  {i,k}=Wrapper.NPGc   (wPhi,wPhit,[],[],wy,initSig,opt);
            end
        end

    case lower('tspAddition')
        filename = [mfilename '.mat']; load(filename);
        fprintf('PET Poisson TV example for TSP\n');

        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        spiral=spiral_m5;
        tfocs=tfocs_200_m5;

        % time cost RMSE
        forSave=[count(:),meanOverK(   fbp,'RMSE'),...
            meanOverK(     pnpg_),...
            meanOverK(  pnpg_n0),...
            meanOverK(   spiral),...
            meanOverK(    tfocs),...
            ];
        save('varyCntPET_TV.data','forSave','-ascii');
         
        % mIdx=6 is also good
        mIdx=4; as=1; k=1;
        fields_={'stepSize','RMSE','time','cost'};
        forSave=addTrace(        pnpg_{mIdx,as,k},     [],fields_); %  1- 4
        forSave=addTrace(     pnpg_n0{mIdx,as,k},forSave,fields_); %  5- 8
        forSave=addTrace(      spiral{mIdx,as,k},forSave,fields_); %  9-12
        forSave=addTrace(   pnpg_nInf{mIdx,as,k},forSave,fields_); % 13-16
        forSave=addTrace(   pnpg_n0m0{mIdx,as,k},forSave,fields_); % 17-20
        forSave=addTrace(       tfocs{mIdx,as,k},forSave,fields_); % 21-24
        forSave=addTrace(      pnpgA0{mIdx,as,k},forSave,fields_); % 25-28
        forSave=addTrace(    pnpgG5A0{mIdx,as,k},forSave,fields_); % 29-32
        forSave=addTrace(    pnpgG5Aq{mIdx,as,k},forSave,fields_); % 33-26
        forSave=addTrace(    pnpgGfA0{mIdx,as,k},forSave,fields_); % 37-40
        forSave=addTrace(    pnpgGfAq{mIdx,as,k},forSave,fields_); % 41-44
        save('cost_itrPET_TV.data','forSave','-ascii');

        nn=128;
        xtrue = read_zubal_emis('nx', nn, 'ny', nn);
        idx=5;
        fprintf('  PNPG: %g%%\n',  pnpg_{idx}.RMSE(end)*100);
        fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
        fprintf(' tfocs: %g%%\n', tfocs{idx}.RMSE(end)*100);
        fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
        img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
        img=showImgMask(  pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_TV_pet.png')
        img=showImgMask( tfocs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'tfocs_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'tfocs_TV_pet.png')
        img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_TV_pet.png')
        img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_TV_pet.png')

        idx=4;
        fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
        fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
        fprintf(' tfocs: %g%%\n', tfocs{idx}.RMSE(end)*100);
        fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
        img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
        img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_TV4_pet.png')
        img=showImgMask( tfocs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'tfocs_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'tfocs_TV4_pet.png')
        img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_TV4_pet.png')
        img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_TV4_pet.png')

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
        t=t+1; forSave(1:length(out.difX),t)=out.difX;
        out=  spiralTV_Long;
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        t=t+1; forSave(1:length(out.difX),t)=out.difX;

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
        fprintf('     FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
        fprintf('    NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
        img=npg{idx}.x; mask=npg{idx}.opt.mask;
        img=showImgMask(   npgTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPGTV_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPGTV_pet.png')
        img=showImgMask(spiralTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRALTV_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRALTV_pet.png')
        img=showImgMask(     fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,     'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),     'FBP_pet.png')
        img=showImgMask(    npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,    'NPGs_pet.eps','psc2'); imwrite(img/max(xtrue(:)),    'NPGs_pet.png')

        idx=4;
        fprintf('   NPGTV: %g%%\n',   npgTV{idx}.RMSE(end)*100);
        fprintf('SPIRALTV: %g%%\n',spiralTV{idx}.RMSE(end)*100);
        fprintf('     FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
        fprintf('    NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
        img=npg{idx}.x; mask=npg{idx}.opt.mask;
        img=showImgMask(   npgTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPGTV_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPGTV_pet2.png')
        img=showImgMask(spiralTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRALTV_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRALTV_pet2.png')
        img=showImgMask(     fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,     'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),     'FBP_pet2.png')
        img=showImgMask(    npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,    'NPGs_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),    'NPGs_pet2.png')

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
        fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
        fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
        fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
        fprintf('NPGscS: (%g%%, %g%%)\n',    npgsc_s.RMSE(end)*100,rmseTruncate(npgsc_s   )*100);
        img=npg{idx}.x; mask=npg{idx}.opt.mask;
        img=showImgMask(   npg{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPG_pet.png')
        img=showImgMask(  npgc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGc_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGc_pet.png')
        img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet.png')
        img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet.png')
        img=showImgMask(  npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGs_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGs_pet.png')
        img=showImgMask( npgsc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGsc_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGsc_pet.png')
        img=showImgMask( npgsc_s.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGscS_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGscS_pet.png')

        idx=4;
        fprintf('   NPG: %g%%\n',   npg{idx}.RMSE(end)*100);
        fprintf('  NPGc: %g%%\n',  npgc{idx}.RMSE(end)*100);
        fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
        fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
        fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
        fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
        img=npg{idx}.x; mask=npg{idx}.opt.mask;
        img=showImgMask(   npg{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPG_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPG_pet2.png')
        img=showImgMask(  npgc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGc_pet2.png')
        img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet2.png')
        img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet2.png')
        img=showImgMask(  npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGs_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGs_pet2.png')
        img=showImgMask( npgsc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGsc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGsc_pet2.png')

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

function [a,b,c]=meanOverK(method,field)
    if(nargin==2)
        a=mean(Cell.getField(method,field),3);
    else
        a=mean(Cell.getField(method,'time'),3);
        b=mean(Cell.getField(method,'cost'),3);
        c=mean(Cell.getField(method,'RMSE'),3);
        a=[a b c];
    end
end
function forSave=addTrace(method,forSave,fields)
    if(~exist('fields','var'))
        fields={'time','cost','RMSE'};
    end
    n=length(fields);
    for i=1:n
        data(:,i)=reshape(getfield(method,fields{i}),[],1);
    end
    forSave=appendColumns(data,forSave);
end
function forSave = appendColumns(col,forSave)
    [r,c]=size(forSave);
    forSave(1:size(col,1),c+1:c+size(col,2))=col;
end



