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

        for k=1:K
            for i=1:length(count)
                [y,Phi,Phit,Psi,Psit,fbpfunc,opt]=loadPET(count(i),opt);

                fbp{i,1,k}.alpha=maskFunc(fbpfunc(y),opt.mask~=0);
                fbp{i,1,k}.RMSE=sqrNorm(fbp{i,1,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);

                fprintf('fbp RMSE=%f\n',sqrNorm(fbp{i,1,k}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha));
                fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
                %u_max=pNorm(Psit(Phit(1-y./opt.bb(:))),inf);
                %u_max=pNorm(Psit(Phit(y)),inf);

                initSig=max(fbp{i,1,k}.alpha,0);

                if(k==2) save(filename); return; end

                opt.fullcont=false;
                j=1;
                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);
                opt.contShrnk=0.1;

                opt.u = 10^a(i)*u_max; opt.proximal='wvltADMM';
                npg   {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgc  {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.proximal='wvltLagrangian';
                spiral{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);

                opt.u = 10^as(i)*u_max; opt.proximal='wvltADMM';
                npgs  {i,j,k}=Wrapper.NPGs   (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgsc {i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);

                opt.u = 10^atv(i)*u_max; opt.proximal='tviso';
                npgTV {i,j,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgTVc{i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                spiralTV{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);

                save(filename);

%               opt.fullcont=true; opt.maxItr=1e4; opt.thresh=1e-6;

%               % for isotv
%               u_max=1;
%               aa =(3:-0.5:-6);
%               opt.u=(10.^aa)*u_max; opt.proximal='tviso';
%               npgTVFull{i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt);
%               for j=1:length(aa); if(aa(j)>-2)
%                   opt.proximal='tviso';
%                   opt.u=10^aa(j)*u_max;
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
%                   opt.proximal='wvltLagrangian';
%                   opt.u=10^aa(j)*u_max;
%                   if(j==1)
%                       spiralFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
%                   else
%                       spiralFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,spiralFull{i,j-1,k}.alpha,opt);
%                   end
%               end; end

%               save(filename);

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

    case 'plot'
        filename = [mfilename '.mat'];
        load(filename);
        fprintf('PET Poisson example\n');

        count = [1e4 1e5 1e6 1e7 1e8 1e9];

        K = 1:3;

        npgTime   = mean(Cell.getField(   npg(:,1,K),'time'),3);
        npgcTime  = mean(Cell.getField(  npgc(:,1,K),'time'),3);
        npgsTime  = mean(Cell.getField(  npgs(:,1,K),'time'),3);
        npgscTime = mean(Cell.getField( npgsc(:,1,K),'time'),3);
        spiralTime= mean(Cell.getField(spiral(:,1,K),'time'),3);

        npgCost   = mean(Cell.getField(   npg(:,1,K),'cost'),3);
        npgcCost  = mean(Cell.getField(  npgc(:,1,K),'cost'),3);
        npgsCost  = mean(Cell.getField(  npgs(:,1,K),'cost'),3);
        npgscCost = mean(Cell.getField( npgsc(:,1,K),'cost'),3);
        spiralCost= mean(Cell.getField(spiral(:,1,K),'cost'),3);

        fbpRMSE   = mean(Cell.getField(   fbp(:,1,K),'RMSE'),3);
        npgRMSE   = mean(Cell.getField(   npg(:,1,K),'RMSE'),3);
        npgcRMSE  = mean(Cell.getField(  npgc(:,1,K),'RMSE'),3);
        npgsRMSE  = mean(Cell.getField(  npgs(:,1,K),'RMSE'),3);
        npgscRMSE = mean(Cell.getField( npgsc(:,1,K),'RMSE'),3);
        spiralRMSE= mean(Cell.getField(spiral(:,1,K),'RMSE'),3);

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
        legend('NPG','NPGs','NPG-TV');


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
        semilogy(m,    npgTime(:,aIdx),'r-*' ); hold on;
        semilogy(m,   npgsTime(:,aIdx),'c-p' );
        semilogy(m, spiralTime(:,aIdx),'k-^' );
        semilogy(m,   npgcTime(:,aIdx),'k*-.');
        semilogy(m,  npgscTime(:,aIdx),'bs-.');
        semilogy(m,spiral8Time(:,aIdx),'go-.');
        legend('npg','npgs','spiral','npgc','npgsc','spiral8');

        forSave=[npgTime(:,aIdx), npgsTime(:,aIdx), npgcTime(:,aIdx), npgscTime(:,aIdx), spiralTime(:,aIdx), spiral8Time(:,aIdx),...
            npgCost(:,aIdx), npgsCost(:,aIdx), npgcCost(:,aIdx), npgscCost(:,aIdx), spiralCost(:,aIdx), spiral8Cost(:,aIdx),...
            npgRMSE(:,aIdx), npgsRMSE(:,aIdx), npgcRMSE(:,aIdx), npgscRMSE(:,aIdx), spiralRMSE(:,aIdx), spiral8RMSE(:,aIdx),...
            m(:)];
        save('varyMeasurementPoisson.data','forSave','-ascii');


end

end
