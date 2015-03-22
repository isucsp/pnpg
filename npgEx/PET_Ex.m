function PET_Ex(op)
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

% PET example, with background noise b
% Vary the total counts of the measurements, with continuation

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        % PET example
        filename = [mfilename '.mat'];
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

%               opt.fullcont=true;
%               opt.u=(10.^aa)*u_max; opt.maxItr=1e4; opt.thresh=1e-12;
%               npgFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgFull{i,k};
%               fprintf('k=%d, good a = 1e%g\n',k,max((aa(out.contRMSE==min(out.contRMSE)))));
%               opt.fullcont=false;

                j=1;
                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example_003',i,j,k);
                opt.u = 10^a(i)*u_max;

                if(i<5) continue;  end

                npg   {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                keyboard
%               spiral{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgc  {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                if(i==5 && k==2)
                    npgsc_s=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);
                end

                opt.u = 10^as(i)*u_max;
%               npgs  {i,j,k}=Wrapper.NPGs   (Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgsc {i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);

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
                save(filename);
            end
        end

    case 'plot'
        filename = [mfilename '.mat'];
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

        forSave=[]; t=0; mIdx=5;
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

        save('cost_itrPET.data','forSave','-ascii');
        mincost=reshape(forSave(:,[1,4,7]),[],1); 
        mincost=min(mincost(mincost~=0));

        figure;
        semilogy(forSave(:,3),forSave(:,1)-mincost,'r'); hold on;
        semilogy(forSave(:,6),forSave(:,4)-mincost,'g');
        semilogy(forSave(:,9),forSave(:,7)-mincost,'b');
        semilogy(forSave(:,12),forSave(:,10)-mincost,'k-.');
        semilogy(forSave(:,15),forSave(:,13)-mincost,'c:');
        legend('npgc','npg','spiral','npgs','npgsc');
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



end
