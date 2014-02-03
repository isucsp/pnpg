%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Mon 03 Feb 2014 12:05:18 AM CST

clear;
setupPath
i=0; NOverP=[]; NOverPi=[];

description='';

subDir='/CircleMask'; %'/cvxHull'; %'/TightMask'; %
imageName='castSim' %'phantom' %'twoMaterials'; %'realct'; %'pellet'; %
runList=[20];
spark=0;

matname='realCTResultForFST_3logspan17points.mat';
%matname='castSimResultForFST_130pointSim3logspan17points.mat';
%matname='castSimForFST_Cut40TightMas130pointSim3logspan17points.mat';

PhiPhitMode='basic'; %'filtered'; %'weighted'; %
saveImg=0;

opt.spectBasis = 'b0';

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
    opt.skipAlpha=1;

    trueAlpha=Img2D(maskIdx);
    trueAlpha=trueAlpha/norm(trueAlpha);
    opt.trueAlpha=trueAlpha;

    opt.maxIeSteps = 1;

    %opt.t3=0;       % set t3 to ignore value of opt.a
    opt.numCall=1;
    aArray=-6.8:0.2:-6.2;
    muLustigArray=logspace(-15,-6,5);
    j=3;
    opt.muLustig=muLustigArray(j);
    opt.skipIe=0;
    %string='realCTCircularMaskFullProjASSparseWithResampleMu1000itrRunOnGridforA_2';
    %string='castSimCircularMaskFullProjASSparseSparkInputWithResampleMu1000';
    %string='twoMaterialsCircularMaskFullProjASWithResampleMu1000itr';
    
    string='castSimCircularMaskFullProj';
    if((~isfield(opt,'t3')) || opt.t3~=0 )
        string=[string 'Sparse'];
    end
    if(spark && ~streq(imageName,'realct'))
        string=[string 'SparkInput'];
    end
    if(opt.numCall>1)
        string=[string 'WithResampleMu' num2str(opt.numCall) 'Times'];
    end
    if(opt.skipIe) string=[string 'KnownIe']; end
    if(opt.maxIeSteps>1) string=[string 'FullAS']; end
    string = [string opt.spectBasis];
    string=[string num2str(opt.maxItr)];
    aArray=-6.5;

    if(opt.numCall>1)
        for i=1:length(aArray)
            opt.a=aArray(i);
            out{i}=beamhardenASSparseResampleMu(B,Bt,C,Ct,y,initSig,opt);
            RMSE(i)=1-(out{i}(end).alpha'*opt.trueAlpha/norm(out{i}(end).alpha))^2;
        end
    else
        for i=1:length(aArray)
            opt.a=aArray(i);
            out{i}=beamhardenSpline(B,Bt,C,Ct,y,initSig,opt);
            RMSE(i)=1-(out{i}.alpha'*opt.trueAlpha/norm(out{i}.alpha))^2;
        end
    end
    if(~ exist('st.mat'))
        save('st.mat','st');
    end
    clear 'st'
    save(string);
    %[out]=beamhardenlog(B,Bt,C,Ct,y,initSig,opt);
   
    % ADD SPARSE RECONSRUCTION 

end
end
%save(matname);

