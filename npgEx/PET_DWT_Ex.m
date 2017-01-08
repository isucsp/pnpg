function PET_DWT_Ex(op)
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
    clear -regexp '(?i)proxopt'
    clear -regexp '(?i)C'
    clear -regexp '(?i)proximal'
    filename = [mfilename '.mat'];
    OPT.mask=[]; OPT.outLevel=1;
    OPT.maxItr=1e3; OPT.thresh=1e-9; OPT.debugLevel=2; OPT.noiseType='poisson';
    OPT.maxItr=1e4; OPT.thresh=1e-9; OPT.debugLevel=1; OPT.noiseType='poisson';
    C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
    PROXOPT.Lip=@(u)u^2; PROXOPT.initStep='fixed';
    PROXOPT.adaptiveStep=false; PROXOPT.backtracking=false;
    proximal=sparseProximal(Psi,Psit,C.prox,'pnpg',PROXOPT);

    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    K=1;

    as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
    a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];
    atv= [-0.5,-0.5,  0,   0, 0.5,   1];

    for k=1:K
    for i=5
        j=1;
        [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
        NLL=@(x) Utils.poissonModel(x,Phi,Phit,y,OPT.bb);

        fbp{i,1,k}.x=fbpfunc(y);
        fbp{i,1,k}.RMSE=sqrNorm(fbp{i,1,k}.x-OPT.trueX)/sqrNorm(OPT.trueX);

        fprintf('fbp RMSE=%f\n',sqrNorm(fbp{i,1,k}.x-OPT.trueX)/sqrNorm(OPT.trueX));
        fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
        u_max=1;
        OPT.u = u_max*10.^a(i);

        initSig=C.prox(fbp{i,1,k}.x);

        fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);

        % BEGIN experiment region,  to delete in the end
        % END experiment region,  to delete in the end

        opt=OPT;
        opt.sigma=10^-5; opt.tau=opt.sigma; opt.maxItr=opt.maxItr*4;
        cpdwt2 {i,j,k}=CP_DWT(Phi,Phit,y,2,Psi,Psit,C,initSig,opt);

        opt=OPT;
        opt.sigma=10^-4; opt.tau=opt.sigma; opt.maxItr=opt.maxItr*4;
        cpdwt1 {i,j,k}=CP_DWT(Phi,Phit,y,1,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.thresh=1e-9;
        opt.L=1/pnpg_{i,j,k}.stepSize(end);
        H.exact=true;
        H.val=@(s) opt.u*norm(s(:),1);
        H.proxConj=@(a,v) max(min(a,opt.u),-opt.u);
        condat   {i,j}=pds(NLL,H,Psi,Psit,C,opt.L,initSig,opt);

        opt=OPT;
        pnpg_   {i,j,k}=pnpg(NLL,proximal,initSig,opt);

        opt=OPT;
        % this version maybe able to used to achieve the smallest objective
        proxOpt=PROXOPT;
        opt.proximal=sparseProximal(Psi,Psit,C.prox,'admm',proxOpt);
        pnpg_admm   {i,j,k}=pnpg(NLL,opt.proximal,initSig,opt);


        opt=OPT; opt.dualGap=true;
        proxOpt=PROXOPT;  proxOpt.dualGap=opt.dualGap;
        opt.proximal=sparseProximal(Psi,Psit,C.prox,'pnpg',proxOpt);
        pnpg_d   {i,j,k}=pnpg(NLL,opt.proximal,initSig,opt);

        opt=OPT; opt.innerThresh=1e-5;
        spiral_m5 {i,j,k}=Wrapper.SPIRAL  (Phi,Phit,Psi,Psit,y,initSig,opt);
        opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-5;
        tfocs_200_m5 {i,j,k}=Wrapper.tfocs    (Phi,Phit,Psi,Psit,y,initSig,opt);

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
        mysave;

        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
        pnpg_n0 {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.adaptiveStep=false;
        pnpg_nInf{i,j,k}=pnpg(NLL,proximal,initSig,opt);

        mysave
    end
    end

case lower('plot')
    filename = [mfilename '.mat']; load(filename);
    fprintf('PET Poisson l1 example\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    K = 5;
          fbp=      fbp(:,:,1:K);
        pnpg_=    pnpg_(:,:,1:K);
       spiral=   spiral_m5(:,:,1:K);
      pnpg_n0=  pnpg_n0(:,:,1:K);
    pnpg_nInf=pnpg_nInf(:,:,1:K);
        pnpgc=    pnpgc(:,:,1:K);
        tfocs=tfocs_200_m5(:,:,1:K);

    figure;
    loglog(count,meanOverK(      fbp,'RMSE'),'b-o'); hold on;
    loglog(count,meanOverK(     pnpg,'RMSE'),'r-*'); hold on;
    loglog(count,meanOverK(   spiral,'RMSE'),'k*-.');
    loglog(count,meanOverK(  pnpg_n0,'RMSE'),'c>-');
    loglog(count,meanOverK(pnpg_nInf,'RMSE'),'gs-');
    loglog(count,meanOverK(    pnpgc,'RMSE'),'bp-.');
    loglog(count,meanOverK(    tfocs,'RMSE'),'kh-');
    legend('fbp','pnpg','spiral','pnpg\_n0','pnpg\_nInf','pnpgc','tfocs');

    figure;
    loglog(count,meanOverK(     pnpg,'time'),'r-*'); hold on;
    loglog(count,meanOverK(   spiral,'time'),'k*-.');
    loglog(count,meanOverK(  pnpg_n0,'time'),'c>-');
    loglog(count,meanOverK(pnpg_nInf,'time'),'gs-');
    loglog(count,meanOverK(    pnpgc,'time'),'bp-.');
    loglog(count,meanOverK(    tfocs,'time'),'kh-');
    legend('pnpg','spiral',' pnpg\_n0','pnpg\_nInf','pnpgc','tfocs');

    % time cost RMSE
    forSave=[count(:),meanOverK(   fbp,'RMSE'),...
        meanOverK(     pnpg),...
        meanOverK(   spiral),...
        meanOverK(  pnpg_n0),...
        meanOverK(pnpg_nInf),...
        meanOverK(    pnpgc),...
        meanOverK(    tfocs),...
        ];
    save('varyCntPET.data','forSave','-ascii');

    keyboard

    % mIdx=6 is also good
    mIdx=5; as=1; k=1;
    fields={'stepSize','RMSE','time','cost'};
    forSave=addTrace(         npg{mIdx,as,k},     [],fields); %  1- 4
    forSave=addTrace(       pnpg_{mIdx,as,k},forSave,fields); %  5- 8
    forSave=addTrace(      spiral{mIdx,as,k},forSave,fields); %  9-12
    forSave=addTrace(   pnpg_nInf{mIdx,as,k},forSave,fields); % 13-16
    forSave=addTrace(     pnpg_n0{mIdx,as,k},forSave,fields); % 17-20
    forSave=addTrace(       tfocs{mIdx,as,k},forSave,fields); % 21-24
    forSave=addTrace(      pnpgA0{mIdx,as,k},forSave,fields); % 25-28
    forSave=addTrace(    pnpgG5A0{mIdx,as,k},forSave,fields); % 29-32
    forSave=addTrace(    pnpgG5Aq{mIdx,as,k},forSave,fields); % 33-26
    forSave=addTrace(    pnpgGfA0{mIdx,as,k},forSave,fields); % 37-40
    forSave=addTrace(    pnpgGfAq{mIdx,as,k},forSave,fields); % 41-44
    save('cost_itrPET.data','forSave','-ascii');
    mincost=reshape(forSave(:,[4,8,12,16,20,24]),[],1); 
    mincost=min(mincost(mincost~=0));

    figure; semilogy(forSave(:,5),'r'); hold on;
    semilogy(forSave(:,13),'b');
    semilogy(forSave(:,17),'k');
    %semilogy(forSave(:,9),'g');
    title('step size versus number of iterations');
    legend('pnpg','npg nInf','pnpg n0');

    figure;
    semilogy(forSave(:, 3),forSave(:, 4)-mincost,'r'); hold on;
    semilogy(forSave(:, 7),forSave(:, 8)-mincost,'g');
    semilogy(forSave(:,11),forSave(:,12)-mincost,'b');
    semilogy(forSave(:,15),forSave(:,16)-mincost,'k');
    semilogy(forSave(:,19),forSave(:,20)-mincost,'c');
    semilogy(forSave(:,23),forSave(:,24)-mincost,'k--');
    legend('npg n4','pnpg n4','spiral','pnpg nInf','pnpg n0','tfocsAT');
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
    semilogy(forSave(:,23),forSave(:,22),'k--');
    legend('npg n4','pnpg n4','spiral','pnpg nInf','pnpg n0','tfocsAT');

    nn=128;
    xtrue = read_zubal_emis('nx', nn, 'ny', nn);
    % attenuation map
    mumap = read_zubal_attn('nx', nn, 'ny', nn);
    imwrite(xtrue/max(xtrue(:)),'pet.png');
    imwrite(mumap/max(mumap(:)),'mumap.png');

    idx=5;
    fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
    fprintf(' PNPGc: %g%%\n', pnpgc{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
    img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
    img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_pet.png')
    img=showImgMask( pnpgc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'PNPGc_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'PNPGc_pet.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet.png')

    idx=4;
    fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
    fprintf(' PNPGc: %g%%\n', pnpgc{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
    img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
    img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_pet2.png')
    img=showImgMask( pnpgc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'PNPGc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)), 'PNPGc_pet2.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet2.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet2.png')


    paperDir='~/research/myPaper/asilomar2014/';
    decide=input(sprintf('start to copy to %s [y/N]?',paperDir),'s');
    if strcmpi(decide,'y')
        system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
    end
    system('rm *_pet.png *_pet.eps *_pet2.eps *_pet2.png');
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



