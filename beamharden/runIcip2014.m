function [conf,opt] = runIcip2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.2 $ $Date: Sun 23 Mar 2014 10:15:25 PM CDT
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0 || ~isempty(runList))
    filename = [mfilename '.mat'];
    if( exist(filename,'file') ) load(filename);
    else save(filename,'filename');
    end
end

if(nargin==0) runList = [0];
elseif(isempty(runList))
    conf=ConfigCT();
    opt = conf.setup();
end

%%%%%%%%%%%%%%%%%%%%%%%%
if(any(runList==0)) % reserved for debug and for the best result
    conf=ConfigCT();
    conf.prjFull = 360/2;
    conf.prjNum = conf.prjFull/2;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.imageName='phantom_1'; %'castSim'; %'phantom' %'twoMaterials'; 
    opt=conf.setup();
    opt.maxIeSteps=2e3;
    opt.maxItr=2000;
    opt.debugLevel=3;
    opt.u = 1e-3;
    initSig = conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    opt.alphaStep = 'FISTA_ADMM_NNL1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    out0=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out0','-append');
end

if(any(runList==0.01)) % reserved for debug and for the best result
    [conf, opt] = defaultInit();
    i=1; j=1;
    conf.prjFull = 80;
    conf.prjNum = conf.prjFull/2;
    conf.imgSize = 1024;
    conf.prjWidth = 1024;
    conf.imageName='castSim'; %'phantom' %'twoMaterials'; 'phantom_1'; %
    %'realct'; 'pellet'; %
    opt.muLustig=logspace(-15,-6,5);
    opt.muLustig=opt.muLustig(3); 3.1623e-11;
    opt.spectBasis = 'dis';
    opt.skipIe=true;
    %opt.continuation = true;
    opt.debugLevel=1;
    opt.maxIeSteps = 100;
    opt=conf.setup(opt);
    initSig = conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    for i=1:3
        opt.u = 10^(-i);
        opt.alphaStep='FISTA_ADMM_NNL1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
        out_knownIe_FISTA{i}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out_knownIe_FISTA','-append');
        opt.alphaStep = 'NCG_PR';
        out_knownIe_NCG_PR{i}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out_knownIe_NCG_PR','-append');
    end
end

if(any(runList==0.1)) % reserved for debug and for the best result
    [conf, opt] = defaultInit();
    i=1; j=1;
    opt.muLustig=3.1623e-11;
    opt.spectBasis = 'b1';
    opt.a=opt.a+log10(0.5);
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=conf.FBP(conf.y); initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    out01=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out01','-append');
end

if(any(runList==0.2)) % reserved for debug and for the best result
    [conf, opt] = defaultInit();
    opt.muLustig=3.1623e-11;
    opt.spectBasis = 'b1';
    opt=conf.setup(opt);
    i=1; j=1;
    fprintf('%s, i=%d, j=%d\n','beamharden',i,j);
    initSig=conf.FBP(conf.y); initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    out02=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out02','-append');
end

if(any(runList==1)) % dis, single AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==2))     % FPCAS
    [conf, opt] = defaultInit();
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
end

if(any(runList==3)) %solve by Back Projection
    [conf, opt] = defaultInit();
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
end

if(any(runList==4))     % FPCAS after linearization
    [conf, opt] = defaultInit();
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
end

if(any(runList==5)) %solve by Back Projection after linearization
    [conf, opt] = defaultInit();
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
end

if(any(runList==6)) % b0, known Ie,
    [conf, opt] = defaultInit();
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
end

if(any(runList==7)) % b1, known Ie,
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    opt.spectBasis = 'b1';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out7{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out7','-append');
    end
end

if(any(runList==71)) % b1, known Ie,
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    conf.PhiMode = 'cpuPrj';
    opt.spectBasis = 'b1';
    opt.skipIe=true;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard known Ie';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out71{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out71','-append');
    end
end

if(any(runList==8)) % reserved for debug and for the best result
    [conf, opt] = defaultInit();
    i=1; j=1;
    opt.spectBasis = 'b1';
    %opt.skipIe = true;
    conf.PhiMode='cpuFanPar'; %'basic'; %'filtered'; %'weighted'; %
    conf.PhiModeGen='cpuFanPar'; %'filtered'; %'weighted'; %
    %opt.maxIeSteps = 100;
    %conf.theta = (0:6:179)';
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    out8{1}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out8','-append');

    conf.PhiMode='basic'; %'filtered'; %'weighted'; %'cpuFanPar'; %
    conf.PhiModeGen='basic'; %'cpuFanPar'; %'filtered'; %'weighted'; %
    opt=conf.setup(opt);
    prefix='BeamHard';
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    %initSig = opt.trueAlpha;
    out8{2}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out8','-append');
end

if(any(runList==11)) % dis, single AS step,
    [conf, opt] = defaultInit();
    intval = 6:-1:1; j=1;
    for i=6:length(intval)
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
end

if(any(runList==12)) % dis, max AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==13)) % b0, single AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==14)) % b0, max AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==15)) % b1, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out15{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out15','-append');
    end
end

if(any(runList==151)) % b1, single AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    conf.PhiMode = 'cpuPrj';
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out151{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out151','-append');
    end
end

if(any(runList==152)) % b1, single AS step,
    [conf, opt] = defaultInit();
    i=1; j=1;
    opt.spectBasis = 'b1';
    opt=conf.setup(opt);
    prefix='BeamHard';
    opt.maxItr = 8000;
    fprintf('%s, i=%d, j=%d\n',prefix,i,j);
    initSig = out15{end}.alpha;
    opt.Ie = out15{end}.Ie;
    out152{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out152','-append');
end

if(any(runList==16)) % b1, max AS step,
    [conf, opt] = defaultInit();
    opt.spectBasis = 'b1';
    opt.maxIeSteps = 100;
    intval = 6:-1:1; j=1;
    for i=1:length(intval)
        conf.prjFull = 360/intval(i);
        conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        prefix='BeamHard';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        initSig=conf.FBP(conf.y);
        initSig = initSig(opt.mask~=0);
        out16{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out16','-append');
    end
end

if(any(runList==21)) % dis, single AS step,
    [conf, opt] = defaultInit();
    intval = 6:-1:1;
    aArray=[-6.5, -9:-4];
    for j=4:length(aArray)
        opt.a = aArray(j);
        for i=1:length(intval)
            conf.prjFull = 360/intval(i);
            conf.prjNum = conf.prjFull/2;
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
end

if(any(runList==22)) % dis, max AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==23)) % b0, single AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==24)) % b0, max AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==25)) % b1, single AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==26)) % b1, max AS step,
    [conf, opt] = defaultInit();
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
end

if(any(runList==30)) % Huber function test for best muHuber
    [conf, opt] = defaultInit();
    j=1;
    opt.spectBasis = 'dis';
    opt=rmfield(opt,'muLustig');
    muHuber=[logspace(-1,-15,8), 0.2];
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    for i=2:length(muHuber)
        opt.muHuber=muHuber(i);
        opt=conf.setup(opt);
        initSig=out30{i-1,j}.alpha;
        opt.Ie=out30{i-1,j}.Ie;
        prefix='Huber function for different mu''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out30{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out30','-append');
    end
end

if(any(runList==31)) % test for different u for lustig
    [conf, opt] = defaultInit();
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
end

if(any(runList==32)) % Huber function test for best muHuber
    [conf, opt] = defaultInit();
    j=1;
    opt.spectBasis = 'dis';
    opt=rmfield(opt,'muLustig');
    muHuber=[1e-16 logspace(-1,-15,8), 1e-16];
    muHuber=muHuber(end:-1:1);
    opt.maxItr=1e3;
    opt=conf.setup(opt);
    initSig=conf.FBP(conf.y);
    initSig = initSig(opt.mask~=0);
    for i=1:length(muHuber)
        opt.muHuber=muHuber(i);
        opt=conf.setup(opt);
        prefix='Huber function for different mu''s';
        fprintf('%s, i=%d, j=%d\n',prefix,i,j);
        out32{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        initSig=out32{i,j}.alpha;
        opt.Ie=out32{i,j}.Ie;
        save(filename,'out32','-append');
    end
end

if(any(runList==999))
    % Test difference between FISTA_Simplex and ActiveSet for IeStep
    conf=ConfigCT();
    conf.prjFull = 360/2;
    conf.prjNum = conf.prjFull/2;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.imageName='phantom_1'; %'castSim'; %'phantom' %'twoMaterials'; 
    opt=conf.setup();
    opt.maxIeSteps=20;
    opt.maxItr=2000;
    opt.debugLevel=5;
    opt.showImg=true;
    initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
    
    opt.alphaStep = 'FISTA_ADMM_NNL1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    for i=3:4
        opt.u=10^(-(i+2));
        out999{i}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out999','-append');
    end
    % Conclusion: Although the two methods converge to the same place, the
    % FISTA version needs too many iterations to reach the point, while
    % Active Set converges quickly. The reason must be the large condition
    % number of the hessian matrix w.r.t. Ie.
end
end

