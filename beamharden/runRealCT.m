%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Fri 24 Jan 2014 11:15:55 AM CST

clear;
i=0; NOverP=[]; NOverPi=[];

description='';

subDir='/CircleMask'; %'/cvxHull'; %'/TightMask'; %
imageName='realct'; %'pellet'; %'castSim' %'phantom' %'twoMaterials'; %
runList=[3,20];
spark=0;

matname='realCTResultForFST_3logspan17points.mat';
%matname='castSimResultForFST_130pointSim3logspan17points.mat';
%matname='castSimForFST_Cut40TightMas130pointSim3logspan17points.mat';

PhiPhitMode='basic'; %'filtered'; %'weighted'; %
saveImg=0;

i=i+1;
configRealCT
NOverP=[NOverP; N/m];
NOverPi=[NOverPi; N/p_I];

%%%%%%%%%%%%%%%%%%%%%%%%
if(sum(runList==3)==1)
    %solve by Back Projection
    j=0;
    prefix='BackProj';
    j=j+1;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    x_BackProj=FBP(y);
    res_BackProj(i,j)=norm(y-Phi(x_BackProj))/normy;
    Img2D_BackProj=reshape(x_BackProj,[my mx]);
%   Img2D_BackProjOrig=Img2D_BackProj;
%   if(size(Img2D_BackProj,1)>size(Img2D,1))
%       temp=(size(Img2D_BackProjOrig,1)-size(Img2D,1))/2;
%       Img2D_BackProj=Img2D_BackProjOrig(temp+1:end-temp,temp+1:end-temp);
%   end
    PSNR_BackProj(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_BackProj(maskIdx),scaleM);
    PSNR_BackProj(i,j,2)=psnr(Img2D,Img2D_BackProj,scale);
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_BackProj');
        imwrite(Img2D_BackProj,...
            [dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.jpg'],'jpg');
    end
    %save([dirname '/' prefix '.mat'],'res_BackProj','PSNR_BackProj');
%   figure; showImg(Img2D_BackProj);
end

if(sum(runList==20)==1)
    j=0;
    prefix='BeamHard';
    %system(['cp ''' which('beamhardenlog') ''' ' dirname '/']);
    x_BackProj=FBP(y);
    res_BackProj=norm(y-Phi(x_BackProj))/normy;
    Img2D_BackProj=reshape(x_BackProj,[my mx]);

    A=HM; At=HMt;
    B=PhiM; Bt=PhiMt;
    C=PsiM; Ct=PsiMt;
    thresh=1e-12;
    maxitr=2e3;
for jj=7 %[5:6] %1:length(rCoeff)
    j=j+1;
    rInit=rCoeff(jj);
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=x_BackProj(:)*1+0*Mask(:)/2; %Img2D; %
    rInit=min(rInit,length(At(y)));
    if(length(At(y))~=m)
        initSig=initSig(maskIdx);
    end
    loadNIST
    massAttenCoeffIon(:,1)=massAttenCoeffIon(:,1)*1e3;
    epsilon=20:150;
    opt.trueMu=interp1(massAttenCoeffIon(10:end,1),...
        massAttenCoeffIon(10:end,2),epsilon,'spline');
    opt.trueIe=gampdf((epsilon-20)*16/100,5,1);
    opt.epsilon=epsilon;
    if(spark)
        opt.trueIe(45)=opt.trueIe(45)*1000;
        opt.trueIe(22)=opt.trueIe(22)*1000;
    end
    opt.trueIe=opt.trueIe/sum(opt.trueIe);
    %opt.mu=trueMu; opt.sampleMode='assigned';
    %opt.muRange=[0.03 30]; opt.sampleMode='exponential';
    opt.logspan=3; opt.sampleMode='logspan';
    opt.K=2; opt.E=17; opt.maxItr=maxitr; opt.rInit=rInit;
    opt.thresh=thresh; opt.mask=Mask;
    opt.useSparse=0; opt.showImg=1; opt.visible=1;
    opt.skipAlpha=0;

    trueAlpha=Img2D(maskIdx);
    trueAlpha=trueAlpha/norm(trueAlpha);
    opt.trueAlpha=trueAlpha;

    runbh
    %[out]=beamhardenlog(B,Bt,C,Ct,y,initSig,opt);
end
end
%save(matname);

% FPC_AS
if(sum(runList==1)==1)
    j=0;
    prefix='FPCAS';
    A=PhiM; At=PhiMt;
    AO=A_operator(A,At);
for a=logspace(-12,-35,10) %[1e-10] %
    j=j+1;
    mu=a*max(abs(At(y)));
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    opt=[];
    load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    opt.x0=reshape(Wt(Img2D_BackProj),m,1); %zeros(m,1); %
    %if(length(At(y))~=m) opt.x0=opt.x0(wvltIdx); end
%   opt.mxitr=maxitr;
%   opt.gtol=thresh*1e4;
%   opt.gtol_scale_x=thresh;
    tic;
    [s_FPCAS, outTmp] = FPC_AS(length(At(y)),AO,y,mu,[],opt);
    t_FPCAS(i,j)=toc;
    obj_FPCAS(i,j)=(norm(y-A(s_FPCAS))^2)/2+mu*sum(abs(s_FPCAS));
    res_FPCAS(i,j)=norm(y-A(s_FPCAS))/normy;
    Img2D_FPCAS=showImgMask(s_FPCAS,Mask); %W(s_FPCAS);
    PSNR_FPCAS(i,j,1)=psnr(Img2D(maskIdx),Img2D_FPCAS(maskIdx),scaleM);
    PSNR_FPCAS(i,j,2)=psnr(Img2D,Img2D_FPCAS,scale);
    a_FPCAS(i,j)=a;
    out_FPCAS{i,j}=outTmp;
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_FPCAS','y');
        imwrite(Img2D_FPCAS,[dirname '/Img2D_' prefix '_' num2str(i) '_'...
            num2str(j) '.jpg'],'jpg');
    end
    save([dirname '/' prefix '.mat'],'a_FPCAS',...
        't_FPCAS','out_FPCAS','res_FPCAS','obj_FPCAS','NOverP','NOverPi',...
        'PSNR_FPCAS');
end
end

% MASK FPC_AS
if(sum(runList==2)==1)
    j=0;
    prefix='maskFPCAS';
    a=[1e-5]; %logspace(-14,-20,7); %
for jj=1:length(a);
    j=j+1;
    A=HM; At=HMt;
    mm_maskFPCAS(i,j)=length(At(y));
    mu=a(jj)*max(abs(At(y)));
    AO=A_operator(A,At);
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    opt=[];
    opt.x0=reshape(Wt(Img2D_BackProj),m,1); %zeros(m,1); %
    if(length(At(y))~=m) opt.x0=opt.x0(wvltIdx); end
%   opt.mxitr=100000;
%   opt.gtol=1e-10;
%   opt.gtol_scale_x=opt.gtol*1e-4;
%   opt.zero=0;
    fprintf('mm=%d, mu=%g\n',mm_maskFPCAS(i,j),mu);
    tic;
    [s_maskFPCAS, outTmp] = FPC_AS(mm_maskFPCAS(i,j),AO,y,mu,[],opt);
    t_maskFPCAS(i,j)=toc;
    obj_maskFPCAS(i,j)=...
        (norm(y-A(s_maskFPCAS))^2)/2+mu*sum(abs(s_maskFPCAS));
    res_maskFPCAS(i,j)=norm(y-A(s_maskFPCAS))/normy;
    outTmp.r=sum(s_maskFPCAS~=0);

    Img2D_maskFPCAS=zeros(my,mx);
    Img2D_maskFPCAS(maskIdx)=PsiM(s_maskFPCAS);

    PSNR_maskFPCAS(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_maskFPCAS(maskIdx),scaleM);
    PSNR_maskFPCAS(i,j,2)=psnr(Img2D,Img2D_maskFPCAS,scale);
    a_maskFPCAS(i,j)=a(jj);
    out_maskFPCAS{i,j}=outTmp;
    if(saveImg)
        save([dirname subDir ...
            '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.mat'], 'Img2D_maskFPCAS','y');
        imwrite(Img2D_maskFPCAS,[dirname subDir ...
            '/Img2D_' prefix '_' num2str(i) '_' num2str(j) '.jpg'],'jpg');
    end
    save([dirname subDir '/' prefix '.mat'],'a_maskFPCAS',...
        'out_maskFPCAS','t_maskFPCAS','res_maskFPCAS','NOverP','NOverPi',...
        'obj_maskFPCAS','mm_maskFPCAS','PSNR_maskFPCAS');
end
end

if(sum(runList==4)==1) %solve by Landweber
    system(['cp ''' which('matrixInv') ''' ' dirname '/']);
    j=0;
    prefix='maskLandweberPhiM';
    j=j+1;
    thresh=1e-8;
    maxItr=3e4;
    z_init=zeros(m,1); %Img2D_LSQR(:); %
%   z_true=Img2D(:);
    A=PhiM; At=PhiMt;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    zLen=length(At(y));
    if(zLen~=m)
        z_init=z_init(maskIdx);
%       z_true=z_true(maskIdx);
    end
    tic;
    [z_MaskBackProj,cntTmp,resTmp,gapTmp]=matrixInv(A,At,y,z_init,...
        thresh,'itrlim',maxItr,'trueVal',[]);
    t_MaskBackProj(i,j)=toc;
    cnt_MaskBackProj(i,j)=cntTmp;
    res_MaskBackProj{i,j}=resTmp;
%   gap_MaskBackProj{i,j}=20*log10(scaleM*sqrt(p_M)./gapTmp);
    if(zLen~=m)
        temp=zeros(m,1);
        temp(maskIdx)=z_MaskBackProj;
        z_MaskBackProj=temp;
    end
    Img2D_MaskBackProj=reshape(z_MaskBackProj,[my mx]);
    PSNR_MaskBackProj(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_MaskBackProj(maskIdx),scaleM);
    PSNR_MaskBackProj(i,j,2)=psnr(Img2D,Img2D_MaskBackProj,scale);
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_MaskBackProj','y');
    end
    save([dirname '/' prefix '.mat'],...  
        't_MaskBackProj','cnt_MaskBackProj','res_MaskBackProj',...
        'NOverP','NOverPi','PSNR_MaskBackProj');
end

if(sum(runList==5)==1)
    j=0;
    system(['cp ''' which('NewDORE') ''' ' dirname '/']);
    %solve by DORE
    prefix='DORE';
    mc=expectHHt(Phi,Phit,N,m,'Phi');
    stepMode=2; %110; %002;
    prefix=sprintf('%s%03d',prefix,stepMode);
    A=H; At=Ht;
    thresh=1e-14;
    maxitr=1e3;
for jj=6:7%:length(rCoeff)
    j=j+1;
    rInit=rCoeff(jj);
    rInit=min(rInit,length(At(y)));
%   rInit=floor(rCoeff(jj)*length(find(abs(wCoeff(:))>0)));
%   load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    initSig=zeros(m,1); %reshape(Wt(Img2D_BackProj),m,1); %
    fprintf('%s i=%d, j=%d\n',prefix,i,j);
    tic;
    [s_DORE,cntTmp,muTmp,outTmp]=...
        NewDORE(A,At,y,rInit,...
        'Thresh',thresh,'visibility',1,'InitialSig',initSig,...
        'maxItr',maxitr,'stepMode',stepMode,'step',1/c);
    t_DORE(i,j)=toc;
    res_DORE(i,j)=norm(y-A(s_DORE))/normy;
    mu_DORE(i,j)=muTmp;
    cnt_DORE(i,j)=cntTmp;
    out_DORE{i,j}=outTmp;
    z_DORE=s_DORE+Ht(y-H(s_DORE))/c;
    s_DORE=reshape(s_DORE,[my mx]);
    Img2D_DOREs=W(s_DORE);
    z_DORE=reshape(z_DORE,[my,mx]);
    Img2D_DOREz=W(z_DORE);

%   Img2D_DOREOrig=Img2D_DOREs;
%   if(size(Img2D_DOREs,1)>size(Img2D,1))
%       temp=(size(Img2D_DOREs,1)-size(Img2D,1))/2;
%       Img2D_DOREs=Img2D_DOREs(temp+1:end-temp,temp+1:end-temp);
%       Img2D_DOREz=Img2D_DOREz(temp+1:end-temp,temp+1:end-temp);
%   end
    
    PSNR_DOREz(i,j,1)=psnr(Img2D(maskIdx),Img2D_DOREz(maskIdx),scaleM);
    PSNR_DOREz(i,j,2)=psnr(Img2D,Img2D_DOREz,scale);
    PSNR_DOREs(i,j,1)=psnr(Img2D(maskIdx),Img2D_DOREs(maskIdx),scaleM);
    PSNR_DOREs(i,j,2)=psnr(Img2D,Img2D_DOREs,scale);
    r_Dore(i,j)=rInit;
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j) ...
            '.mat'],'Img2D_DOREz','Img2D_DOREs','y');
        imwrite(Img2D_DOREs,...
            [dirname '/Img2D_' prefix 's_' num2str(i) '_' num2str(j)...
            '.jpg'],'jpg');
    end
    save([dirname '/' prefix '.mat'],'out_DORE','mu_DORE',...
        'cnt_DORE','t_DORE','res_DORE','r_Dore','NOverP',...
        'NOverPi','PSNR_DOREs','PSNR_DOREz');
end
end

if(sum(runList==6)==1)
    %solve by maskDORE
    j=0;
    system(['cp ''' which('NewDORE') ''' ' dirname '/']);
    prefix='maskDORE';
    stepMode=112; %112; %
    prefix=sprintf('%s%03d',prefix,stepMode);
    if(mod(stepMode,10)==2)
        mc=2; %expectHHt(PhiM,PhiMt,N,p_M,'PhiM');
    elseif(mod(stepMode,10)==3)
        mc=repmat(diagonal(:),Num_proj,1);
    else mc=1; end
%   load([dirname '/' prefix '.mat']);
%   rInit=floor(length(find(abs(wCoeff(:))>0)));
    A=HM; At=HMt;
    thresh=1e-14;
    maxitr=1e3;
for jj=2:3 %[5:6] %1:length(rCoeff)
    j=j+1;
    rInit=rCoeff(jj);
    fprintf('%s i=%d, j=%d\n',prefix,i,j);
    %load([dirname '/Img2D_BackProj_1_1.mat']);
    initSig=zeros(m,1); %reshape(Wt(Img2D_BackProj),m,1); %
    rInit=min(rInit,length(At(y)));
    if(length(At(y))~=m)
        initSig=initSig(wvltIdx);
    end
    [s_maskDORE,cntTmp,muTmp,outTmp]=...
        NewDORE(A,At,y,rInit,...
        'Thresh',thresh,'visibility',1,'InitialSig',initSig,...
        'maxItr',maxitr,'stepMode',stepMode,'step',1./mc);
    t_maskDORE(i,j)=toc;
    cnt_maskDORE(i,j)=cntTmp;
    mu_maskDORE{i,j}=muTmp;
    out_maskDORE{i,j}=outTmp;
    resTmp=y-A(s_maskDORE);  res_maskDORE(i,j)=norm(resTmp)/normy;
    if(isvector(muTmp)) z_maskDORE=s_maskDORE+At(muTmp.*resTmp);
    else z_maskDORE=s_maskDORE+At(muTmp*resTmp); end

    Img2D_maskDOREs=zeros(my,mx);
    Img2D_maskDOREs(maskIdx)=PsiM(s_maskDORE);
    Img2D_maskDOREz=zeros(my,mx);
    Img2D_maskDOREz(maskIdx)=PsiM(z_maskDORE);

    PSNR_maskDOREs(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_maskDOREs(maskIdx),scaleM);
    PSNR_maskDOREs(i,j,2)=psnr(Img2D,Img2D_maskDOREs,scale);
    PSNR_maskDOREz(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_maskDOREz(maskIdx),scaleM);
    PSNR_maskDOREz(i,j,2)=psnr(Img2D,Img2D_maskDOREz,scale);
    r_maskDore(i,j)=rInit;
    if(saveImg)
        save([dirname subDir ...
            '/Img2D_' prefix '_' num2str(i) '_' num2str(j) '.mat'],...
            'Img2D_maskDOREz','Img2D_maskDOREs','y');
        imwrite(Img2D_maskDOREs,[dirname subDir ...
            '/Img2D_' prefix 's_' num2str(i) '_' num2str(j) '.jpg'],'jpg');
    end
    if(0)
        afun=@(xxx,yyy) A_LSQR(xxx,PhiM,PhiMt,yyy);
        [w,flgTmp,relresTmp,cntTmp,resTmp]=lsqr(afun,resTmp,thresh,maxitr,...
            [],[],[]);
        temp=zeros(m,1);
        temp(maskIdx)=w;
        w=temp;
        Img2D_maskDOREzeb=Img2D_maskDOREs+reshape(w,[my mx]);
        PSNR_maskDOREzeb(1,i,j)=...
            psnr(Img2D(maskIdx),Img2D_maskDOREzeb(maskIdx),scaleM);
        PSNR_maskDOREzeb(2,i,j)=psnr(Img2D,Img2D_maskDOREzeb,scale);
        save([dirname subDir ...
            '/Img2D_' prefix '_' num2str(i) '_' num2str(j) '.mat'],...
            'Img2D_maskDOREzeb','PSNR_maskDOREzeb','-append');
    end
    save([dirname subDir '/' prefix '.mat'],...
        'cnt_maskDORE','res_maskDORE','t_maskDORE','r_maskDore',...
        'mu_maskDORE','NOverP','NOverPi','out_maskDORE',...
        'PSNR_maskDOREs','PSNR_maskDOREz');
end
end

if(sum(runList==9)==1)
    %solve by Mask Iterative Non-negative Threshold
    j=1;
    tic;
    thresh=1e-14;
    A=PhiM; At=PhiMt;
    load([dirname '/Img2D_BackProj_1_1.mat']);
    xInit=Img2D_BackProj(:); %zeros(m,1); %
    if(length(At(y))~=m)
        xInit=xInit(maskIdx);
    end
    [Img2D_maskINT,Count_maskINT]=intDore(A,At,y,xInit,thresh);
    t_maskINT(i,j)=toc;
    temp=zeros(size(Mask));
    temp=temp(:);
    temp(maskIdx)=Img2D_maskINT;
    Img2D_maskINT=reshape(temp,[my mx]);
    PSNR_maskINT(i,j,1)=psnr(Img2D(maskIdx),Img2D_maskINT(maskIdx),scaleM);
    PSNR_maskINT(i,j,2)=psnr(Img2D,Img2D_maskINT,scale);
    save([dirname '/maskINT.mat'],'PSNR_maskINT',...
                                                        't_maskINT');
    if(saveImg)
        save([dirname '/Img2D_maskINT_' num2str(i) '_' num2str(j)...
                                            '.mat'],'Img2D_maskINT');
    end
end

if(sum(runList==11)==1)
    j=0;
    prefix='LSQRPhiM0';
    afun=@(xxx,yyy) A_LSQR(xxx,PhiM,PhiMt,yyy);
    %load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    initSig=zeros(m,1); %reshape(Wt(Img2D_LSQR),m,1); %Img2D_BackProj(:); %
    thresh=1e-14;
for maxitr=[2e4]
    j=j+1;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    if(length(afun(y,'transp'))~=m)
        initSig=initSig(maskIdx);
    end
    tic;
    % initSig here in lsqr function is useless
    [x_LSQR,flgTmp,relresTmp,cntTmp,resTmp]=lsqr(afun,y,thresh,maxitr,...
        [],[],initSig);
    t_LSQR(i,j)=toc;
    relres_LSQR(i,j)=relresTmp;
    cnt_LSQR(i,j)=cntTmp;
    flg_LSQR(i,j)=flgTmp;
    res_LSQR{i,j}=resTmp/normy;
    temp=zeros(m,1);
    temp(maskIdx)=x_LSQR;
    x_LSQR=temp;
    Img2D_LSQR=reshape(x_LSQR,[my mx]);
    PSNR_LSQR(i,j,1)=psnr(Img2D(maskIdx),Img2D_LSQR(maskIdx),scaleM);
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j) ...
            '.mat'],'Img2D_LSQR','y');
    end
    save([dirname '/' prefix '.mat'],'t_LSQR','flg_LSQR',...
        'res_LSQR','cnt_LSQR','relres_LSQR','NOverP','NOverPi',...
        'PSNR_LSQR');
end
end

if(sum(runList==12)==1)
    %solve by Mask Back Projection
    tic;
    s_init=zeros(length(HMt(y)),1);
    [s_MaskBackProj,Count]=matrixInv(HM,HMt,y,s_init,thresh);
    t_MaskBackProj(i,j)=toc;
    temp=zeros(m,1);
    temp(wvltIdx)=s_MaskBackProj;
    s_MaskBackProj=reshape(temp,my, mx);
    Img2D_MaskBackProj=W(s_MaskBackProj).*Mask;
    PSNR_MaskBackProj(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_MaskBackProj(maskIdx),scaleM);
    save([dirname '/MaskBackProj_' num2str(i) '_' num2str(j) '.mat'],...
        'PSNR_MaskBackProj','t_MaskBackProj',...
        'Count','NOverP','NOverPi','resTmp','PSNR_MaskBackProj');
    if(saveImg)
        save([dirname '/Img2D_MaskBackProj_' num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_MaskBackProj','y');
    end
end

if((sum(runList==13)==1)) %solve by GPSR_BB
    j=0;
    prefix='GPSR';
    A=H; At=Ht;
%   reshape(Wt(Img2D_BackProj),m,1);
    tol=1e-5;
    load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    initSig=reshape(Wt(Img2D_BackProj),m,1); %zeros(m,1); %
    if(length(At(y))~=m) initSig=initSig(wvltIdx); end
%for a=1e-4; %logspace(-15,-5,11)
    a=logspace(-12,-3,10) %[1e-7] %
for jj=8;
    j=j+1;
    fprintf('%s i=%d, j=%d\n',prefix,i,j);
    tau=a(jj)*max(abs(At(y)));
    tic;
    [s_GPSR,s_GPSRdb,objTmp,times,debTmp,mesTmp]=...
        GPSR_BB(y,A,tau,'AT',At,'Debias',1,...
        'ToleranceA',tol,...
        'Verbose',0,'Initialization',initSig);
    t_GPSR(i,j)=toc;
    if(debTmp~=0) res_GPSR(i,j)=norm(y-A(s_GPSRdb))/normy;
    else res_GPSR(i,j)=norm(y-A(s_GPSR))/normy; end
    obj_GPSR{i,j}=objTmp;
    deb_GPSR(i,j)=debTmp;
    mes_GPSR{i,j}=mesTmp;
    a_GPSR(i,j)=a(jj);
    if isempty(s_GPSRdb) s_GPSRdb=zeros(size(s_GPSR)); end
    if(length(s_GPSR)~=m)
        temp=zeros(m,1);
        temp(wvltIdx)=s_GPSR;
        s_GPSR=temp;
        temp=zeros(m,1);
        temp(wvltIdx)=s_GPSRdb;
        s_GPSRdb=temp;
    end
    s_GPSR=reshape(s_GPSR,my,mx);
    s_GPSRdb=reshape(s_GPSRdb,my,mx);
    Img2D_GPSR=W(s_GPSR);
    Img2D_GPSRdb=W(s_GPSRdb);
    PSNR_GPSR(i,j,1)=psnr(Img2D(maskIdx),Img2D_GPSR(maskIdx),scaleM);
    PSNR_GPSR(i,j,2)=psnr(Img2D,Img2D_GPSR,scale);
    PSNR_GPSRdb(i,j,1)=psnr(Img2D(maskIdx),Img2D_GPSRdb(maskIdx),scaleM);
    PSNR_GPSRdb(i,j,2)=psnr(Img2D,Img2D_GPSRdb,scale);
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_GPSRdb','Img2D_GPSR','y');
    end
    save([dirname '/' prefix '.mat'],...
        'obj_GPSR','t_GPSR','deb_GPSR','mes_GPSR',...
        'NOverP','NOverPi', 'res_GPSR','a_GPSR',...
        'PSNR_GPSRdb','PSNR_GPSR');
end
end

if((sum(runList==14)==1)) %solve by GPSR_BB
    j=0;
    prefix='maskGPSR';
    A=HM; At=HMt;
    tol=1e-5;
%for a=logspace(-12,-3,10)
    a=logspace(-12,-3,10); %[1e-10] %
for jj=8 %length(a);
    j=j+1;
    fprintf('%s i=%d, j=%d\n',prefix,i,j);
    load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    initSig=reshape(Wt(Img2D_BackProj),m,1); %zeros(m,1); %
    if(length(At(y))~=m) initSig=initSig(wvltIdx); end
    tic;
    tau_maskGPSR(i,j)=a(jj)*max(abs(At(y)));
    [s_maskGPSR,s_maskGPSRdb,objTmp,times,debTmp,mesTmp]=...
        GPSR_BB(y,A,tau_maskGPSR(i,j),'AT',At,'Debias',1,...
        'ToleranceA',tol,...
        'Verbose',0,'Initialization',initSig);
    t_maskGPSR(i,j)=toc;
    obj_maskGPSR(i,j)=...
        (norm(y-A(s_maskGPSR))^2)/2+tau_maskGPSR(i,j)*sum(abs(s_maskGPSR));
    if(debTmp~=0) res_maskGPSR(i,j)=norm(y-A(s_maskGPSRdb))/normy;
    else res_maskGPSR(i,j)=norm(y-A(s_maskGPSR))/normy; end
    vecObj_maskGPSR{i,j}=objTmp;
    deb_maskGPSR(i,j)=debTmp;
    mes_maskGPSR{i,j}=mesTmp;
    a_maskGPSR(i,j)=a(jj);

    Img2D_maskGPSR=zeros(my,mx);
    if(length(s_maskGPSR)~=m)
        Img2D_maskGPSR(maskIdx)=PsiM(s_maskGPSR);
        Img2D_maskGPSRdb=Img2D_maskGPSR;
        if(debTmp~=0)
            Img2D_maskGPSRdb(maskIdx)=PsiM(s_maskGPSRdb);
        end
    end

    PSNR_maskGPSR(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_maskGPSR(maskIdx),scaleM);
    PSNR_maskGPSR(i,j,2)=psnr(Img2D,Img2D_maskGPSR,scale);
    PSNR_maskGPSRdb(i,j,1)=...
        psnr(Img2D(maskIdx),Img2D_maskGPSRdb(maskIdx),scaleM);
    PSNR_maskGPSRdb(i,j,2)=psnr(Img2D,Img2D_maskGPSRdb,scale);
    if(saveImg)
        save([dirname subDir ...
            '/Img2D_' prefix '_' num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_maskGPSRdb','Img2D_maskGPSR','y');
    end
    save([dirname subDir '/' prefix '.mat'],...
        'obj_maskGPSR','t_maskGPSR','deb_maskGPSR','mes_maskGPSR',...
        'NOverP','NOverPi', 'res_maskGPSR','a_maskGPSR',...
        'PSNR_maskGPSRdb','PSNR_maskGPSR','vecObj_maskGPSR',...
        'tau_maskGPSR');
end
end

% %Standard Filtered Backprojection
% tic;
% clc;
% display('Now FBP...');
% Img2D_InvRadon=iradon(CTdata_mtrx(:,end:-1:1),theta/pi*180+360-theta(end)/pi*180,'nearest',filter,1,N_freq); %Uniform
% t_InvRadon=toc;
% Img2D_InvRadon_max=max(Img2D_InvRadon(:));

if(sum(runList==15)==1)
    j=0;
    prefix='LSQRPhi';
    afun=@(xxx,yyy) A_LSQR(xxx,Phi,Phit,yyy);
    load([dirname '/Img2D_BackProj_' num2str(i) '_1.mat']);
    initSig=Img2D_BackProj(:); %zeros(m,1); %reshape(Wt(Img2D_LSQR),m,1); %
    thresh=1e-14;
for maxitr=[2e4];
    j=j+1;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    if(length(afun(y,'transp'))~=m)
        initSig=initSig(maskIdx);
    end
    tic;
    % initSig here in lsqr function is useless
    [x_LSQR,flgTmp,relresTmp,cntTmp,resTmp]=lsqr(afun,y,thresh,maxitr,...
        [],[],initSig);
    t_LSQR(i,j)=toc;
    relres_LSQR(i,j)=relresTmp;
    cnt_LSQR(i,j)=cntTmp;
    flg_LSQR(i,j)=flgTmp;
    res_LSQR{i,j}=resTmp/normy;
    Img2D_LSQR=reshape(x_LSQR,[my mx]);
    PSNR_LSQR(i,j,1)=psnr(Img2D(maskIdx),Img2D_LSQR(maskIdx),scaleM);
    if(saveImg)
        save([dirname '/Img2D_' prefix '_' num2str(i) '_' num2str(j) ...
            '.mat'],'Img2D_LSQR','y');
    end
    save([dirname '/' prefix '.mat'],'t_LSQR','flg_LSQR',...
        'res_LSQR','cnt_LSQR','relres_LSQR','NOverP','NOverPi',...
        'PSNR_LSQR');
end
end

if(sum(runList==16)==1)
    j=0;
    system(['cp ''' which('NonnegDORE') ''' ' dirname '/']);
    prefix='NonnegDOREPhi';
for rInit=1;
    j=j+1;
    A=Phi; At=Phit;
    thresh=1e-14;
    maxitr=1e4;
    stepMode=110;
    initSig=zeros(m,1); %reshape(Wt(Img2D_LSQR),m,1);
    fprintf('%s i=%d, j=%d\n',prefix,i,j);
    tic;
    if(length(At(y))~=m)
        initSig=initSig(maskIdx);
    end
    [s_maskDORE,cntTmp,resVecTmp,idxTmp,muTmp,DOREisGoodTmp]=...
        NonnegDORE(A,At,y,rInit,...
        'Thresh',thresh,'visibility',1,'InitialSig',initSig,...
        'maxItr',maxitr,'stepMode',stepMode);
    t_maskDORE(i,j)=toc;
    cnt_maskDORE(i,j)=cntTmp;
    idx_maskDORE{i,j}=idxTmp;
    resVec_maskDORE{i,j}=resVecTmp;
    mu_maskDORE(i,j)=muTmp;
    DOREisGood_maskDORE{i,j}=DOREisGoodTmp;
    resTmp=y-A(s_maskDORE);  res_maskDORE(i,j)=norm(resTmp)/normy;
    z_maskDORE=s_maskDORE+At(resTmp);
    if(length(At(y))~=m)
        temp=zeros(m,1);
        temp(maskIdx)=s_maskDORE;
        s_maskDORE=temp;
        temp(maskIdx)=z_maskDORE;
        z_maskDORE=temp;
    end
    Img2D_maskDOREs=reshape(s_maskDORE,[my mx]);
    Img2D_maskDOREz=reshape(z_maskDORE,[my,mx]);
    r_maskDore(i,j)=rInit;
    if(saveImg)
        save([dirname '/Img2D_' prefix '_'  num2str(i) '_' num2str(j)...
            '.mat'],'Img2D_maskDOREz','Img2D_maskDOREs','y');
    end
    save([dirname '/' prefix '.mat'],'idx_maskDORE','resVec_maskDORE',...
        'cnt_maskDORE','res_maskDORE','t_maskDORE','r_maskDore',...
        'mu_maskDORE','NOverP','NOverPi','DOREisGood_maskDORE');
end
end

%solve by ECME
if(any(runList==17))
    display('Now ECME...');
    system(['cp ''' which('ECMEls') ''' ' dirname '/']);
    prefix=['ECME'];
    rInit=7000;
    A=HM; At=HMt;
    thresh=1e-14;
    %load([dirname '/Img2D_BackProj_1_1.mat']);
    initSig=zeros(m,1); %reshape(Wt(Img2D_BackProj),m,1); %
    rInit=min(rInit,length(At(y)));
    if(length(At(y))~=m)
        initSig=initSig(wvltIdx);
    end
    mu_ECME=0.0023192;
    load('parPrj/diagnal512.mat');
    mu_ECME=1./(diagnal(:)+250);
    temp=[zeros(1,180); diagnal1];
    mu_ECME=[diagnal(:)+0 temp(:)];
    tic;
    [s_ECME,cnt_ECME,out_ECME]=ECMEls(A,At,y,initSig,rInit,thresh,mu_ECME);
    t_ECME=toc; r_ECME=rInit;
    if(size(mu_ECME,2)==1)
        if(isscalar(mu))
            z_ECME=s_ECME+At(y-A(s_ECME))/mu_ECME;
        else
            z_ECME=s_ECME+At((y-A(s_ECME))./mu_ECME);
        end
    else
        temp=y-A(s_ECME);
        % CAUTION: temp can not be reused since being changed.
        temp=solveTridiag(mu_ECME(:,2),mu_ECME(:,1),...
            [mu_ECME(2:end,2);0],temp,'inv');
        clear('temp');
        z_ECME=s_ECME+At(temp);
    end
    Img2D_ECME=zeros(my,mx);
    Img2D_ECMEz=zeros(my,mx);
    if(length(At(y))~=m)
        Img2D_ECME(maskIdx)=PsiM(s_ECME);
        Img2D_ECMEz(maskIdx)=PsiM(z_ECME);
    else
        Img2D_ECME(:)=Psi(s_ECME);
        Img2D_ECMEz(:)=Psi(z_ECME);
    end
    if(saveImg)
        save([dirname subDir ...
            '/Img2D_' prefix '.mat'],'Img2D_ECMEz','Img2D_ECME');
    end
    save([dirname subDir '/' prefix '.mat'],...
        'cnt_ECME','out_ECME','t_ECME','r_ECME',...
        'mu_ECME');
end

%solve by ADOREls
if(any(runList==18))
    display('Now ADORE...');
    system(['cp ''' which('ADORE') ''' ' dirname '/']);
    A=HM; At=HMt;
    stepMode=002;
    prefix=['ADORE'];
    Len_thresh=5000;
    thresh=1e-8;
    tic;
    [s_ADORE,out_ADORE,cnt_ADORE,Golden_Iter]=...
        ADORE(A,At,y,'Thresh',thresh,'SearchLen',Len_thresh,...
        'stepMode',stepMode,'step',1/mc);
    t_ADORE=toc;
    z_ADORE=s_ADORE+At(y-A(s_ADORE))/mc;

    Img2D_ADORE=zeros(my,mx);
    Img2D_ADORE(maskIdx)=PsiM(s_ADORE);
    Img2D_ADOREz=zeros(my,mx);
    Img2D_ADOREz(maskIdx)=PsiM(z_ADORE);

    if(saveImg)
        save([dirname subDir '/Img2D_' prefix '.mat'],...
            'Img2D_ADOREz','Img2D_ADORE');
    end
    save([dirname subDir '/' prefix '.mat'],...
        'cnt_ADORE','out_ADORE','Golden_Iter','t_ADORE');
end

%solve by NIHTls
if(any(runList==19))
    j=0;
    display('Now NIHT...');
    system(['cp ''' which('hard_l0_Mterm_kq') ''' ' dirname '/']);
    prefix=['NIHT'];
    A=HM; At=HMt;

    thresh=1e-3;
    maxitr=1e4;
for jj=[1 3]%:length(rCoeff)
    j=j+1;
    r_init=rCoeff(jj);
    tic;
    [s_NIHT,outTmp]=...
        hard_l0_Mterm_kq(y,A,p_I,r_init,thresh,'P_trans',At,'maxIter',maxitr);
    t_NIHT(i,j)=toc;
    r_NIHT(i,j)=r_init;
    cnt_NIHT(i,j)=length(outTmp.resVec);
    resTmp=y-A(s_NIHT);
    res_NIHT(i,j)=norm(resTmp)/normy;
    out_NIHT{i,j}=outTmp;
    z_NIHT=s_NIHT+At(resTmp)/mc;

    Img2D_NIHT=zeros(my,mx);
    Img2D_NIHT(maskIdx)=PsiM(s_NIHT);
    Img2D_NIHTz=zeros(my,mx);
    Img2D_NIHTz(maskIdx)=PsiM(z_NIHT);

    if(saveImg)
        save([dirname subDir '/Img2D_' prefix '_' num2str(i)...
            '_' num2str(j) '.mat'],...
            'Img2D_NIHT','Img2D_NIHTz');
    end
    save([dirname subDir '/' prefix '.mat'],...
        'r_NIHT','res_NIHT','out_NIHT','cnt_NIHT','t_NIHT');
end
end

