function PET_l1_Ex(op)
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
        OPT.maxItr=1e4; OPT.thresh=1e-6; OPT.debugLevel=1; OPT.noiseType='poisson';
        OPT.contShrnk=0.1; OPT.contGamma=15; OPT.noBackground=true;

        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        K=1;

        as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
        a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];

        OPT.mask=[];
        for k=1:K
            for i=length(count):-1:1
                j=1;
                [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);

                fbp{i,1,k}.alpha=maskFunc(fbpfunc(y),OPT.mask~=0);
                fbp{i,1,k}.RMSE=sqrNorm(fbp{i,1,k}.alpha-OPT.trueAlpha)/sqrNorm(OPT.trueAlpha);

                fprintf('fbp RMSE=%f\n',sqrNorm(fbp{i,1,k}.alpha-OPT.trueAlpha)/sqrNorm(OPT.trueAlpha));
                fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
                u_max=1;

                initSig=max(fbp{i,1,k}.alpha,0);

                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);

                opt=OPT; opt.u = 10^a(i)*u_max; opt.proximal='wvltADMM'; opt.adaptiveStep=false;
                pnpg_nInf   {i,j,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u = 10^a(i)*u_max; opt.proximal='wvltADMM'; opt.cumuTol=0; opt.incCumuTol=false;
                pnpg_n0   {i,j,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);

                save(filename);
                continue;

                if any(i==[3])
                    opt=OPT; opt.u = 10^a(i)*u_max; opt.thresh=1e-13;
                    spiral_long{i,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
                    save(filename);
                    continue;

                    opt=OPT; opt.u = 10^a(i)*u_max; opt.thresh=1e-10; opt.proximal='wvltADMM';
                    pnpg_n4_long{i,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npg_n4_long{i,k}=Wrapper.NPG    (Phi,Phit,Psi,Psit,y,initSig,opt);

                    opt=OPT; opt.u = 10^a(i)*u_max; opt.thresh=1e-10; opt.proximal='wvltADMM'; opt.adaptiveStep=false;
                    pnpg_nInf_long{i,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    
                    opt=OPT; opt.u = 10^a(i)*u_max; opt.thresh=1e-10; opt.proximal='wvltADMM'; opt.cumuTol=0; opt.incCumuTol=false;
                    pnpg_n0_long{i,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                else
                    continue;
                end
                save(filename);
                continue;

                opt=OPT; opt.u = 10^a(i)*u_max; opt.proximal='wvltFADMM';
                fpnpg  {i,j,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                opt=OPT; opt.u = 10^a(i)*u_max; opt.proximal='wvltADMM';
                pnpg   {i,j,k}=Wrapper.PNPG    (Phi,Phit,Psi,Psit,y,initSig,opt);
                pnpgc  {i,j,k}=Wrapper.PNPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                npg    {i,j,k}=Wrapper.NPG     (Phi,Phit,Psi,Psit,y,initSig,opt);
                opt=OPT; opt.u = 10^a(i)*u_max;
                spiral {i,j,k}=Wrapper.SPIRAL  (Phi,Phit,Psi,Psit,y,initSig,opt);

                % for wavelet l1 norm
                aa = (3:-0.5:-6);
                opt=OPT; opt.fullcont=true; opt.u=(10.^aa)*u_max; opt.proximal='wvltADMM';
                pnpgFull {i,k}=Wrapper.PNPG (Phi,Phit,Psi,Psit,y,initSig,opt);
                opt=OPT; opt.fullcont=true; opt.u=(10.^aa)*u_max; opt.proximal='wvltFADMM';
                fpnpgFull{i,k}=Wrapper.PNPG (Phi,Phit,Psi,Psit,y,initSig,opt);
                for j=1:length(aa); if(aa(j)>-2)
                    opt=OPT; opt.u=10^aa(j)*u_max; opt.proximal='wvltLagrangian';
                    if(j==1)
                        spiralFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,initSig,opt);
                    else
                        spiralFull{i,j,k}=Wrapper.SPIRAL (Phi,Phit,Psi,Psit,y,spiralFull{i,j-1,k}.alpha,opt);
                    end
                end; end

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

    case lower('plot')
        filename = [mfilename '.mat']; load(filename);
        fprintf('PET Poisson l1 example\n');

        count = [1e4 1e5 1e6 1e7 1e8 1e9];

        K = 1;
          pnpg=  pnpg(:,:,1:K);
           fbp=   fbp(:,:,1:K);
        spiral=spiral(:,:,1:K);
           npg=   npg(:,:,1:K);
         fpnpg= fpnpg(:,:,1:K);
         pnpgc= pnpgc(:,:,1:K);

        figure;
        loglog(count,meanOverK(  pnpg,'RMSE'),'r-*'); hold on;
        loglog(count,meanOverK(   fbp,'RMSE'),'b-o');
        loglog(count,meanOverK(spiral,'RMSE'),'k-^');
        loglog(count,meanOverK(   npg,'RMSE'),'k*-.');
        loglog(count,meanOverK( fpnpg,'RMSE'),'c>-');
        loglog(count,meanOverK( pnpgc,'RMSE'),'gs-');
        legend('pnpg','fbp','spiral','npg','fpnpg','pnpgc');

        figure;
        loglog(count,meanOverK(  pnpg,'time'),'r-*'); hold on;
        loglog(count,meanOverK(spiral,'time'),'k*-.');
        loglog(count,meanOverK(   npg,'time'),'c>-');
        loglog(count,meanOverK( fpnpg,'time'),'gs-');
        loglog(count,meanOverK( pnpgc,'time'),'gs-');
        legend('pnpg','spiral',' npg','fpnpg','pnpgc');

        forSave=[...
            meanOverK(  pnpg),...
            meanOverK(spiral),...
            meanOverK(  pnpg),...
            meanOverK( fpnpg),...
            meanOverK(   fbp,'RMSE'),count(:)];
        save('varyCntPET.data','forSave','-ascii');

        keyboard

        mIdx=3; as=1; k=1;
        fields={'stepSize','RMSE','time','cost'};
        forSave=addTrace(   npg_n4_long{mIdx,as,k},     [],fields); %  1- 4
        forSave=addTrace(  pnpg_n4_long{mIdx,as,k},forSave,fields); %  5- 8
        forSave=addTrace(   spiral_long{mIdx,as,k},forSave,fields); %  9-12
        forSave=addTrace(pnpg_nInf_long{mIdx,as,k},forSave,fields); % 13-16
        forSave=addTrace(  pnpg_n0_long{mIdx,as,k},forSave,fields); % 17-20

        save('cost_itrPET.data','forSave','-ascii');
        mincost=reshape(forSave(:,[4,8,12,16,20]),[],1); 
        mincost=min(mincost(mincost~=0));

        figure; semilogy(forSave(:,5),'r'); hold on;
        semilogy(forSave(:,13),'b');
        semilogy(forSave(:,17),'k');
        %semilogy(forSave(:,9),'g');
        title('step size versus number of iterations');
        legend('pnpg','npg nInf','pnpg n0','spiral');

        figure;
        semilogy(forSave(:, 3),forSave(:, 4)-mincost,'r'); hold on;
        semilogy(forSave(:, 7),forSave(:, 8)-mincost,'g');
        semilogy(forSave(:,11),forSave(:,12)-mincost,'b');
        semilogy(forSave(:,15),forSave(:,16)-mincost,'k');
        semilogy(forSave(:,19),forSave(:,20)-mincost,'c');
        legend('npg n4','pnpg n4','spiral','pnpg nInf','pnpg n0');
        hold on;

        keyboard

        idx=min(find(forSave(:,10)<1e-6));
        plot(forSave(idx,9),forSave(idx,7)-mincost,'bo');
        xxx=idx;
        idx=min(find(forSave(10:end,14)<1e-6))+10;
        plot(forSave(idx,13),forSave(idx,11)-mincost,'k*');
        xxx=[xxx;idx];  xxx=xxx(:)';
        save('cost_itrPETidx.data','xxx','-ascii');

        figure;
        semilogy(forSave(:, 3),forSave(:, 2),'r'); hold on;
        semilogy(forSave(:, 7),forSave(:, 6),'g');
        semilogy(forSave(:,11),forSave(:,10),'b');
        semilogy(forSave(:,15),forSave(:,14),'k');
        semilogy(forSave(:,19),forSave(:,18),'c');
        legend('npg n4','pnpg n4','spiral','pnpg nInf','pnpg n0');

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


        paperDir='~/research/myPaper/asilomar2015/';
        decide=input(sprintf('start to copy to %s [y/N]?',paperDir));
        if strcmpi(decide,'y')
            system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
        end
        system('rm *_pet.eps *_pet2.eps *_pet2.png');
        close all;

    case 'fullplot'
        filename = [mfilename '.mat'];
        load(filename);

        k=1;
        aa =(3:-0.5:-6);
        for i=1:length(count)
            pnpgContRMSE  {i,k} = [  pnpgFull{i,k}.contRMSE(:);  pnpgFull{i,k}.RMSE(end)]; out=pnpgContRMSE{i,k};
            fprintf('i=%d, good a = 1e%g PNPG\n',i,max((aa(out==min(out)))));
            fpnpgContRMSE {i,k} = [ fpnpgFull{i,k}.contRMSE(:); fpnpgFull{i,k}.RMSE(end)]; out=fpnpgContRMSE{i,k};
            fprintf('i=%d, good a = 1e%g FPNPG\n',i,max((aa(out==min(out)))));
            spiralContRMSE {i,k} = Cell.getField(spiralFull(i,:,k),'RMSE'); out=fpnpgContRMSE{i,k};
            fprintf('i=%d, good a = 1e%g SPIRAL\n',i,max((aa(out==min(out)))));
        end

        for i=1:length(count)
            figure;
            semilogy(aa(1:length(pnpgContRMSE{i})),pnpgContRMSE{i},'r-*'); hold on;
            semilogy(aa(1:length(fpnpgContRMSE{i})),fpnpgContRMSE{i},'g-o');
            semilogy(aa(1:length(spiralContRMSE{i})),spiralContRMSE{i},'b-s');
            title(num2str(i));
            legend('PNPG','FPNPG','SPIRAL');
            aaa(i)=min(pnpgContRMSE{i});
            bbb(i)=min(fpnpgContRMSE{i});
            ccc(i)=min(spiralContRMSE{i});
        end
        figure; semilogy(aaa,'r-*'); hold on;
        semilogy(bbb,'g-o');
        semilogy(ccc,'b-s');
        title('rmse vs count');
        legend('PNPG','FPNPG','SPIRAL');
end

end

function [a,b,c]=meanOverK(method,field)
    if(nargin==2)
        a=mean(Cell.getField(method,field),3);
    else
        a=mean(Cell.getField(method,'time'),3);
        b=mean(Cell.getField(method,'cost'),3);
        c=mean(Cell.getField(method,'RMSE'),3);
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



