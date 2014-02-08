function runIcip2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.2 $ $Date: Sat 08 Feb 2014 02:54:00 AM CST
%   v_0.2:      Changed to class oriented for easy configuration

filename = [mfilename '.mat'];
if(~exist(filename,'file')) save(filename,''); end
if(nargin==0) runList=3; end

[conf, opt] = defaultInit();

%%%%%%%%%%%%%%%%%%%%%%%%
if(any(runList==3)) %solve by Back Projection
    j=1; i=1;
    prefix='BackProj';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    opt=conf.setup(opt);
    out3.x_BackProj=conf.FBP(conf.y);
    out3.res_BackProj=norm(conf.y-Phi(x_BackProj))/normy;
    out3.Img2D_BackProj=reshape(x_BackProj,[my mx]);
    out3.PSNR_BackProj(i,j,1)=psnr(Img2D(maskIdx),Img2D_BackProj(maskIdx),scaleM);
    out3.PSNR_BackProj(i,j,2)=psnr(Img2D,Img2D_BackProj,scale);
    save(filename,'out3','-append');
end

if(any(runList==21)) % dis, single AS step,
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out21{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out21','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==22)) % dis, max AS step,
    opt.maxIeSteps = 100;
    intval = 6:-1:1;
    i=1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out22{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out22','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==23)) % b0, single AS step,
    opt.spectBasis = 'b0';
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out23{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out23','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==24)) % b0, max AS step,
    opt.spectBasis = 'b0';
    opt.maxIeSteps = 100;
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out24{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out24','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==25)) % b1, single AS step,
    opt.spectBasis = 'b1';
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out25{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out25','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==26)) % b1, max AS step,
    opt.spectBasis = 'b1';
    opt.maxIeSteps = 100;
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=1:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);
            prefix='BeamHard';
            fprintf('%s, i=%d, j=%d\n',prefix,i,j);
            initSig=conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            out26{j,i}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out26','-append');
        end
    end
    [conf, opt] = defaultInit();
end

% ADD SPARSE RECONSRUCTION 

if(any(runList==20)) % beamhardening with refinement
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
        save(filename,'out20','-append');
    end
end

if(any(runList==0))
    servers={'linux-3', 'research-4'}
    for i=1:length(servers)
        eval(['!scp ' servers{i} ':/local/renliang/imgRecSrc/beamharden/' filename ' temp.mat']);
        st = load('temp.mat');
        names = fieldnames(st);
        for i=1:length(names)
            fprintf('saving %s...\n',names{i});
            eval([names{i} '= st.' names{i} ';']);
            save(filename,names{i});
        end
    end
end

end

% Don't change the default values arbitrarily, it will cause the old code 
% unrunnable
function [conf, opt] = defaultInit()
    conf = ConfigCT();
    conf.maskType='CircleMask'; %'cvxHull'; %'TightMask'; %
    conf.imageName='castSim'; %'phantom' %'twoMaterials'; %'realct'; %'pellet'; %
    conf.PhiPhitMode='basic'; %'filtered'; %'weighted'; %
    conf.spark=0;

    % the higher, the more information. Set to 0 to turn off.
    opt.debugLevel = 7;

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
    opt.showImg=0;
    opt.visible=1;
    opt.skipAlpha=0;
    opt.maxIeSteps = 1;
    %opt.t3=0;       % set t3 to ignore value of opt.a
    opt.numCall=1;
    opt.muLustig=1e-13; % logspace(-15,-6,5);
    opt.skipIe=0;
    opt.a=-6.5;  % aArray=-6.8:0.2:-6.2;
end
