%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Thu 06 Feb 2014 09:57:57 PM CST

clear;
setupPath
filename = [mfilename '.mat'];

conf = ConfigCT();
conf.maskType='CircleMask'; %'cvxHull'; %'TightMask'; %
conf.imageName='castSim'; %'phantom' %'twoMaterials'; %'realct'; %'pellet'; %
conf.PhiPhitMode='basic'; %'filtered'; %'weighted'; %
conf.spark=0;

opt.rInit=5000;
opt.spectBasis = 'dis';
opt.saveImg=0;
opt.thresh=1e-12;
opt.maxItr=2e3;
%opt.muRange=[0.03 30];
opt.logspan=3;
opt.sampleMode='logspan'; %'assigned'; %'exponential'; %
opt.K=2;
opt.E=17;
opt.useSparse=0;
opt.showImg=1;
opt.visible=1;
opt.skipAlpha=0;
opt.maxIeSteps = 100;
%opt.t3=0;       % set t3 to ignore value of opt.a
opt.numCall=1;
opt.muLustig=1e-13; % logspace(-15,-6,5);
opt.skipIe=0;
opt.a=-6.5;  % aArray=-6.8:0.2:-6.2;

runList=[20];

%%%%%%%%%%%%%%%%%%%%%%%%
if(any(runList==3)) %solve by Back Projection
    j=0;
    prefix='BackProj';
    j=j+1;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    opt=conf.setup(opt);
    out3.x_BackProj=conf.FBP(conf.y);
    out3.res_BackProj=norm(conf.y-Phi(x_BackProj))/normy;
    out3.Img2D_BackProj=reshape(x_BackProj,[my mx]);
    out3.PSNR_BackProj(i,j,1)=psnr(Img2D(maskIdx),Img2D_BackProj(maskIdx),scaleM);
    out3.PSNR_BackProj(i,j,2)=psnr(Img2D,Img2D_BackProj,scale);
    save(filename,'out3','-append');
end

if(any(runList==20))
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    out20{i}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    out20{i}.RMSE=1-(out20{i}.alpha'*opt.trueAlpha/norm(out20{i}.alpha))^2;
    save(filename,'out20','-append');
end


% ADD SPARSE RECONSRUCTION 


if(any(runList==21)) % beamhardening with refinement
    j=0;
    prefix='BeamHard';
    x_BackProj=conf.FBP(conf.y);
    for jj=7 %[5:6] %1:length(rCoeff)
        j=j+1;
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=x_BackProj(:)*1+0*Mask(:)/2; %Img2D; %

        aArray=-6.8:0.2:-6.2;
        muLustigArray=logspace(-15,-6,5);
        j=3;
        opt.muLustig=muLustigArray(j);

        aArray=-6.5;
        for i=1:length(aArray)
            opt.a=aArray(i);
            out{i}=beamhardenASSparseResampleMu(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            RMSE(i)=1-(out{i}(end).alpha'*opt.trueAlpha/norm(out{i}(end).alpha))^2;
        end
        save(filename,'out21','-append');
    end
end
