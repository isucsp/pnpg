function [conf,opt] = runIcip2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.2 $ $Date: Thu 27 Mar 2014 10:21:08 PM CDT
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0 || ~isempty(runList))
    filename = [mfilename '.mat'];
    if( exist(filename,'file') ) load(filename);
    else save(filename,'filename');
    end
end

if(nargin==0) runList = [0];
elseif(isempty(runList))
    conf=ConfigCT(); opt = conf.setup(); return;
end

% runList rules, abc
% a:
%       0: castSim
%       9: phantom_1

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

% dis, known Ie,
if(any(runList==001))
    conf=ConfigCT();
    opt.skipIe=true;
    opt.maxItr=1e3;
    opt.showImg=true; opt.debugLevel=5;
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u = 10.^[-1 -2 -3 -4 -5 -6 -7];
    for i=1:1
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        for j=3:7
            fprintf('%s, i=%d, j=%d\n','CPLS',i,j);
            opt.u=u(j);
            out001{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'out001','-append');
            initSig=out001{i,j}.alpha;
        end
    end
end

if(any(runList==002))     % FPCAS
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    aArray=[-10:-4];
    for i=2:3
        for j=1:length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            opt.a = aArray(j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
            opt=conf.setup(opt);

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            mu=10^opt.a*max(abs(At(conf.y)));
            option.x0=conf.FBP(conf.y);
            option.x0 = conf.Psit(option.x0(opt.mask~=0));
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            out002{i,j}=out; out002{i,j}.alpha = conf.Psi(s);
            out002{i,j}.opt = opt;
            out002{i,j}.RMSE=1-(out002{i,j}.alpha'*opt.trueAlpha/norm(out002{i,j}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out002','-append');
        end
    end
end

if(any(runList==003)) %solve by Back Projection
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    for i=1:length(prjFull)
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup();
        out003{i}.img=conf.FBP(conf.y);
        out003{i}.alpha=out003{i}.img(opt.mask~=0);
        out003{i}.RSE=norm(conf.y-conf.Phi(out003{i}.alpha))/norm(conf.y);
        out003{i}.RMSE=1-(out003{i}.alpha'*opt.trueAlpha/norm(out003{i}.alpha)/norm(opt.trueAlpha))^2;
        save(filename,'out003','-append');
    end
end

if(any(runList==004))     % FPCAS after linearization
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    aArray=[-10:-4];
    for i=2:3
        for j=1:length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            opt.a = aArray(j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
            opt=conf.setup(opt);

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
            mu=10^opt.a*max(abs(At(y)));
            option.x0=conf.FBP(y);
            option.x0 = conf.Psit(option.x0(opt.mask~=0));
            [s, out] = FPC_AS(length(At(y)),AO,y,mu,[],option);
            out004{i,j}=out; out004{i,j}.alpha = conf.Psi(s);
            out004{i,j}.opt = opt;
            out004{i,j}.RMSE=1-(out004{i,j}.alpha'*opt.trueAlpha/norm(out004{i,j}.alpha)/norm(opt.trueAlpha))^2;
            save(filename,'out004','-append');
        end
    end
end

if(any(runList==005)) %solve by Back Projection after linearization
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    for i=1:length(prjFull)
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup();
        y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
        out005{i}.img=conf.FBP(y);
        out005{i}.alpha=out005{i}.img(opt.mask~=0);
        out005{i}.RSE=norm(y-conf.Phi(out005{i}.alpha))/norm(y);
        out005{i}.RMSE=1-(out005{i}.alpha'*opt.trueAlpha/norm(out005{i}.alpha)/norm(opt.trueAlpha))^2;
        save(filename,'out005','-append');
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

if(any(runList==009))     % FPCAS after linearization
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u=10.^[-1 -2 -3 -4 -5 -6 -7];
    for i=1:6
        for j=1:7
            fprintf('%s, i=%d, j=%d\n','SPIRAL-TAP',i,j);
            conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
            opt.u = u(j);
            opt=conf.setup(opt);
            opt.maxItr=2e3;
            opt.thresh=1e-12;

            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection

            initSig = conf.FBP(conf.y);
            initSig = initSig(opt.mask~=0);
            subtolerance=1e-6;
            out=[];
            [out.alpha, out.p, out.cost, out.reconerror, out.time, ...
                out.solutionpath] = ...
                SPIRALTAP_mod(conf.y,conf.Phi,opt.u,'penalty','ONB',...
                'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
                'initialization',initSig,'maxiter',opt.maxItr,...
                'miniter',0,'stopcriterion',3,...
                'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                'subtolerance',subtolerance,'monotone',1,...
                'saveobjective',1,'savereconerror',1,'savecputime',1,...
                'reconerrortype',2,...
                'savesolutionpath',1,'verbose',100);
            out.opt=opt; out009{i,j}=out;
            save(filename,'out009','-append');
        end
    end
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

% dis, max 20 (default) AS steps, CastSim, FISTA_ADMM_NNL1(default)
if(any(runList==012))
    conf=ConfigCT();
    opt.debugLevel=5; opt.showImg=true;
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    for i=1:1
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=1e-6;
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

% Test difference between FISTA_Simplex and ActiveSet for IeStep
if(any(runList==999))
    conf=ConfigCT();
    conf.prjFull = 360/2;
    conf.prjNum = conf.prjFull/2;
    conf.imgSize = 256;
    conf.prjWidth = 256;
    conf.imageName='phantom_1'; %'castSim'; %'phantom' %'twoMaterials'; 
    opt=conf.setup();
    opt.maxIeSteps=20;
    opt.maxItr=2000;
    opt.debugLevel=1;
    opt.showImg=true;
    initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
    
    opt.alphaStep = 'FISTA_ADMM_NNL1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
    opt.uMode='abs';
    % 1:3, where i=2 (i.e. u=1e-3 gives the best solution
    for i=1:-3
        opt.u=10^(-(i+1));
        out999{i}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out999','-append');
    end
    opt.uMode='relative';
    % 4:6, where i=4 (i.e. opt.a=1e-6, -> u=1.58e-3) gives the best 
    % solution. However, this round has opt.u not updated while running.
    % 7, opt=1e-6, with updating opt.u, not as good as i=4;
    opt.continuation=true;
    opt.uMode='abs';
    opt.maxIeSteps=1;
    opt.skipIe=true;
    for i=9:9
        opt.u=2e-3;
        out999{i}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out999','-append');
    end
    save(filename,'out999','-append');

    % Conclusion: Although the two methods converge to the same place, the
    % FISTA version needs too many iterations to reach the point, while
    % Active Set converges quickly. The reason must be the large condition
    % number of the hessian matrix w.r.t. Ie.
end
end

