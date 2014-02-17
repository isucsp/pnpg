function runIcip2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.2 $ $Date: Mon 17 Feb 2014 01:23:06 AM CST
%   v_0.2:      Changed to class oriented for easy configuration

filename = [mfilename '.mat'];
if(~exist(filename,'file'))
    save(filename,'filename');
else
    load(filename);
end
if(nargin==0) runList=0; end

[conf, opt] = defaultInit();

%%%%%%%%%%%%%%%%%%%%%%%%
if(any(runList==0)) % reserved for debug and for the best result
    i=1; j=1;
    opt.spectBasis = 'b1';
    opt.skipIe = true;
    %opt.maxIeSteps = 100;
    %conf.theta = (0:6:179)';
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    %initSig=conf.FBP(conf.y);
    %initSig = initSig(opt.mask~=0);
    initSig = opt.trueAlpha;
    keyboard
    out0=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out0','-append');
    [conf, opt] = defaultInit();
end

if(any(runList==1)) % dis, single AS step,
    intval = 6:-1:1; j=1;
    opt.spectBasis = 'dis';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out1{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out1','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==2))     % FPCAS
    intval = 6:-1:1;
    aArray=[-10:-4];
    for i=1:length(intval)
        for j=1:length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            opt.a = aArray(j);
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            mu=10^opt.a*max(abs(At(conf.y)));
            option.x0=conf.FBP(conf.y);
            option.x0 = conf.Psit(option.x0(opt.mask~=0));
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            out2{i,j}=out; out2{i,j}.alpha = conf.Psi(s);
            out2{i,j}.opt = opt;
            out2{i,j}.RMSE=1-(out2{i,j}.alpha'*opt.trueAlpha/norm(out2{i,j}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out2','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==3)) %solve by Back Projection
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        out3{i}.img=conf.FBP(conf.y);
        out3{i}.alpha=out3{i}.img(opt.mask~=0);
        out3{i}.RSE=norm(conf.y-conf.Phi(out3{i}.alpha))/norm(conf.y);
        out3{i}.RMSE=1-(out3{i}.alpha'*opt.trueAlpha/norm(out3{i}.alpha)/norm(opt.trueAlpha))^2;
        save(filename,'out3','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==4))     % FPCAS after linearization
    intval = 6:-1:1;
    aArray=[-10:-4];
    for i=1:length(intval)
        for j=1:length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            opt.a = aArray(j);
            conf.theta = (0:intval(i):179)';
            opt=conf.setup(opt);

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
            mu=10^opt.a*max(abs(At(y)));
            option.x0=conf.FBP(y);
            option.x0 = conf.Psit(option.x0(opt.mask~=0));
            [s, out] = FPC_AS(length(At(y)),AO,y,mu,[],option);
            out4{i,j}=out; out4{i,j}.alpha = conf.Psi(s);
            out4{i,j}.opt = opt;
            out4{i,j}.RMSE=1-(out4{i,j}.alpha'*opt.trueAlpha/norm(out4{i,j}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out4','-append');
        end
    end
    [conf, opt] = defaultInit();
end

if(any(runList==5)) %solve by Back Projection after linearization
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
        out5{i}.img=conf.FBP(y);
        out5{i}.alpha=out5{i}.img(opt.mask~=0);
        out5{i}.RSE=norm(y-conf.Phi(out5{i}.alpha))/norm(y);
        out5{i}.RMSE=1-(out5{i}.alpha'*opt.trueAlpha/norm(out5{i}.alpha)/norm(opt.trueAlpha))^2;
        save(filename,'out5','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==6)) % b0, known Ie,
    intval = 6:-1:1; j=1;
    opt.spectBasis = 'b0';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out6{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out6','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==7)) % b1, known Ie,
    intval = 6:-1:1; j=1;
    opt.spectBasis = 'b1';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out7{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out7','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==11)) % dis, single AS step,
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out11{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out11','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==12)) % dis, max AS step,
    opt.maxIeSteps = 100;
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out12{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out12','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==13)) % b0, single AS step,
    opt.spectBasis = 'b0';
    intval = 6:-1:1;
    j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out13{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out13','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==14)) % b0, max AS step,
    opt.spectBasis = 'b0';
    opt.maxIeSteps = 100;
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out14{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out14','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==15)) % b1, single AS step,
    opt.spectBasis = 'b1';
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out15{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out15','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==16)) % b1, max AS step,
    opt.spectBasis = 'b1';
    opt.maxIeSteps = 100;
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.theta = (0:intval(i):179)';
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out16{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out16','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==21)) % dis, single AS step,
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=4:length(aArray)
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
        if(j==1) aaa=5; else aaa=1; end
        for i=aaa:length(intval)
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
    for j=2:length(aArray)
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
        if(j==1) aaa=4
        else aaa=1; end
        for i=aaa:length(intval)
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
        if(j==1) aaa=2; else aaa=1; end
        for i=aaa:length(intval)
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
        if(j==1) aaa=3; else aaa=1; end
        for i=aaa:length(intval)
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

if(any(runList==30)) % Huber function test for best muHuber
    j=1;
    opt.spectBasis = 'dis';
    opt=rmfield(opt,'muLustig');
    muHuber=[logspace(-1,-15,8), 0.2];
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    for i=1:length(muHuber)
        opt.muHuber=muHuber(i);
        opt=conf.setup(opt);
        prefix='Huber function for different mu''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out30{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        initSig=out30{i,j}.alpha;
        opt.Ie=out30{i,j}.Ie;
        save(filename,'out30','-append');
    end
    [conf, opt] = defaultInit();
end

if(any(runList==31)) % test for different u for lustig
    j=1;
    opt.spectBasis = 'dis';
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    aArray=[-6.5, -9:-4];
    for i=1:length(aArray)
        opt.a = aArray(i);
        opt=conf.setup(opt);
        prefix='continuous method for different u''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out31{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        initSig=out31{i,j}.alpha;
        opt.Ie=out31{i,j}.Ie;
        save(filename,'out31','-append');
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

if(any(runList==1000))
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
    conf.PhiMode='basic'; %'filtered'; %'weighted'; %
    conf.PhiModeGen='cpuPar'; %'filtered'; %'weighted'; %
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
    opt.a=-6.5;  % aArray=-6.8:0.2:-6.2;
    opt.K=2;
    opt.E=17;
    opt.useSparse=0;
    opt.showImg=0;
    opt.visible=1;
    %opt.t3=0;       % set t3 to ignore value of opt.a
    opt.numCall=1;
    opt.muLustig=1e-13; % logspace(-15,-6,5);
    opt.skipAlpha=0;
    opt.maxAlphaSteps = 1;
    opt.skipIe=0;
    opt.maxIeSteps = 1;
end
