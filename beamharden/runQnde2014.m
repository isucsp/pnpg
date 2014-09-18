function [conf,opt] = runQnde2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polychromatic Sparse Image Reconstruction and Mass Attenuation Spectrum 
%            Estimation via B-Spline Basis Function Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:  Changed to class oriented for easy configuration

if(nargin==0) runList = [0];
elseif(isempty(runList))
    conf=ConfigCT(); opt = conf.setup(); return;
end
paperDir = '~/research/myPaper/qnde2014/doc/';

%%%%%%%%%%%%%%%%%%%%%%%%

% dis discretized, polychromatic model
if(any(runList==001))
    filename = [mfilename '_001.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    conf=ConfigCT('castSim','CircleMask','gpuPrj');
    prjFull = [60, 80, 100, 120, 180, 360];
    for i=5:length(prjFull)
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt.maxItr=2e3; opt.thresh=1e-6; opt.errorType=0;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

        % unknown ι(κ), NPG-AS
        u  =  10.^[-4  -4   -4   -4   -4   -4];
        for j=[3]
            fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
            opt.u=u(i)*10^(j-2);
            npgas{i,j}=BH.NPG_AS(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npgas','-append');
        end
        continue;

        % known ι(κ), 

        j=1;
        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        fbp{i}.img=conf.FBP(-log(conf.y));
        fbp{i}.alpha=fbp{i}.img(opt.mask~=0);
        fbp{i}.RSE=sqrNorm(-log(conf.y)-conf.Phi(fbp{i,j}.alpha))/sqrNorm(conf.y);
        fbp{i}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);
        save(filename,'fbp','-append');

        aArray=[-10:-2];
        for j=5:7
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);


            mu=10^aArray(j)*max(abs(At(conf.y)));

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            option.x0 = conf.Psit(maskFunc(conf.FBP(conf.y),opt.mask~=0));
            option.mxitr=opt.maxItr;
            [s, out] = FPC_AS(length(At(conf.y)),AO,conf.y,mu,[],option);
            fpcas{i,j}=out; fpcas{i,j}.alpha = conf.Psi(s);
            fpcas{i,j}.fVal(1)=0.5*sqrNorm(conf.Phi(fpcas{i,j}.alpha)-conf.y);
            fpcas{i,j}.fVal(2)=sqrNorm(fpcas{i,j}.alpha.*(fpcas{i,j}.alpha<0));
            fpcas{i,j}.fVal(3)=pNorm(conf.Psit(fpcas{i,j}.alpha),1);
            fpcas{i,j}.opt = opt;
            out004{i,j}.RMSE=1-(out004{i,j}.alpha'*opt.trueAlpha)^2/sqrNorm(out004{i,j}.alpha)/sqrNorm(opt.trueAlpha);
            fprintf('fpcas RMSE=%g\n',fpcas{i}.RMSE);
            save(filename,'fpcas','-append');
        end
    end
end

% noiselessly linearized
if(any(runList==002))
    filename = [mfilename '_002.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('opt');
    conf=ConfigCT('castSim','CircleMask','gpuPrj');
    prjFull = [60, 80, 100, 120, 180, 360];
    for i=1:length(prjFull)
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt.maxItr=2e3; opt.thresh=1e-6; opt.errorType=0;
        opt=conf.setup(opt);
        conf.y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);

        u=10.^[-7 -7 -7 -6 -5 -5];
        for j=4:4
            fprintf('%s, i=%d, j=%d\n','linearized NPG',i,j);
            opt.u = u(i)*10^(j-2);
            opt.continuation = false; opt.alphaStep='FISTA_ADMM_NNL1';
            npg{i,j}=BH.LinNPG(conf.Phi,conf.Phit,...
                conf.Psi,conf.Psit,conf.y,initSig,opt);
            save(filename,'npg','-append');
        end
        continue

        fprintf('%s, i=%d, j=%d\n','BackProj',i,j);
        fbp{i}.img=conf.FBP(y);
        fbp{i}.alpha=fbp{i}.img(opt.mask~=0);
        fbp{i}.RSE=sqrNorm(y-conf.Phi(fbp{i,j}.alpha))/sqrNorm(y);
        fbp{i}.RMSE=sqrNorm(fbp{i,j}.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
        fprintf('fbp RMSE=%g\n',fbp{i,j}.RMSE);
        save(filename,'fbp','-append');
       
        aArray=[-10:-2];
        for j=3:6; %length(aArray)
            fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);
            mu=10^aArray(j)*max(abs(At(y)));

            A = @(xx) conf.Phi(conf.Psi(xx));
            At = @(yy) conf.Psit(conf.Phit(yy));
            AO=A_operator(A,At);
            y = conf.Phi(opt.trueAlpha); % equivalent to linear projection


            option.x0 = conf.Psit(maskFunc(conf.FBP(y),opt.mask~=0));
            option.mxitr=opt.maxItr;
            [s, out] = FPC_AS(length(At(y)),AO,y,mu,[],option);
            fpcas{i,j}=out; fpcas{i,j}.alpha = conf.Psi(s);
            fpcas{i,j}.fVal(1)=0.5*sqrNorm(conf.Phi(fpcas{i,j}.alpha)-y);
            fpcas{i,j}.fVal(2)=sqrNorm(fpcas{i,j}.alpha.*(fpcas{i,j}.alpha<0));
            fpcas{i,j}.fVal(3)=pNorm(conf.Psit(fpcas{i,j}.alpha),1);
            fpcas{i,j}.opt = opt;
            out004{i,j}.RMSE=1-(out004{i,j}.alpha'*opt.trueAlpha)^2/sqrNorm(out004{i,j}.alpha)/sqrNorm(opt.trueAlpha);
            fprintf('fpcas RMSE=%g\n',fpcas{i}.RMSE);
            save(filename,'fpcas','-append');
        end

        u=10.^[-6 -5 -5 -5 -5 -5]; opt.u = u(i); j=1;
        fprintf('%s, i=%d, j=%d\n','SPIRAL-TAP',i,j);
        y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
        initSig = maskFunc(conf.FBP(y),opt.mask~=0);
        out=[]; subtolerance=1e-5;
        [out.alpha, out.p, out.cost, out.reconerror, out.time, out.difAlpha] = ...
            SPIRALTAP_mod(y,conf.Phi,opt.u,'penalty','ONB',...
            'AT',conf.Phit,'W',conf.Psi,'WT',conf.Psit,'noisetype','gaussian',...
            'initialization',initSig,'maxiter',opt.maxItr,...
            'miniter',0,'stopcriterion',3,...
            'tolerance',opt.thresh,'truth',opt.trueAlpha,...
            'subtolerance',subtolerance,'monotone',1,...
            'saveobjective',1,'savereconerror',1,'savecputime',1,...
            'reconerrortype',2,'savedifalpha',1,...
            'savesolutionpath',0,'verbose',10);
        out.opt=opt; spiral{i,j}=out;
        spiral{i,j}.fVal(1)=0.5*sqrNorm(conf.Phi(spiral{i,j}.alpha)-y);
        spiral{i,j}.fVal(2)=sqrNorm(spiral{i,j}.alpha.*(spiral{i,j}.alpha<0));
        spiral{i,j}.fVal(3)=sum(abs(conf.Psit(spiral{i,j}.alpha)));
        save(filename,'spiral','-append');
    end
end

% dis, max 20 (default) AS steps, CastSim, FISTA_ADMM_NNL1(default)
% Also compare the continuation and non-continuation
% This section has very strong connection with 001
if(any(runList==012))
    load(filename,'out012');
    conf=ConfigCT();
    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    u  =  10.^[-5  -4   -4   -4   -4   -4];
    opt.continuation=false; opt.maxIeSteps=1;

    opt.saveAnimate=true;
    for i=2:2
        conf.PhiMode='gpuPrj'; opt.thresh=1e-6;
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end
    return;

    j=2; opt.continuation=true;
    for i=1:-6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end

    j=3; opt.continuation=true; opt.skipIe=true;
    for i=2:-6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
    end
    j=4; opt.continuation=false; opt.skipIe=true;
    for i=1:6
        conf.prjFull = prjFull(i); conf.prjNum = conf.prjFull/2;
        opt=conf.setup(opt);
        initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
        opt.u=u(i);
        out012{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
            conf.Psi,conf.Psit,conf.y,initSig,opt);
        save(filename,'out012','-append');
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

if(any(runList==033))
    filename = [mfilename '_033.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    conf.prjFull = 360; conf.prjNum = conf.prjFull/2; opt.u = 1e-4;
    opt=conf.setup(opt); initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
    opt.maxIeSteps=1; opt.thresh=1e-12;
    opt.maxItr=2000;

    opt.skipIe=true;
    opt.alphaStep='FISTA_ADMM_NNL1';
    newnewRestart=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'newnewRestart','-append');
    return;

    opt.skipIe=true;
    i=1; j=1; opt.alphaStep='FISTA_ADMM_NNL1';
    newRestart{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'newRestart','-append');
    return;

    opt.skipIe=true;
    i=2; j=1; opt.alphaStep='FISTA_ADMM_NNL1';
    oldRestart{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'oldRestart','-append');
end
   

% dis, compare the FISTA_ADMM_NNL1 and NCG_PR for both continuation and 
% non-continuation
if(any(runList==021))
    filename = [mfilename '_021.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    conf.prjFull = 360; conf.prjNum = conf.prjFull/2; opt.u = 1e-4;
    opt=conf.setup(opt); initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
    opt.maxIeSteps=1; opt.thresh=1e-6;
    opt.maxItr=2000;
   
    % opt.adaptiveStep=false;
    % temp=opt.thresh; opt.thresh=1e-15;
    % i=1; j=3; opt.alphaStep='IST_ADMM_NNL1';
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % opt.thresh=temp;
    % save(filename,'out','-append');

    % i=1; j=1; opt.alphaStep='FISTA_ADMM_NNL1';
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');


    % opt.debugLevel=6;
    % i=1; j=6; opt.alphaStep='FISTA_ADMM_NNL1';
    % opt.CenterB=false;
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');
    % rmfield(opt,'CenterB');

    % i=1; j=5; opt.alphaStep='FISTA_ADMM_NNL1';
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');

    % % i=1; j=2; opt.alphaStep='NCG_PR';
    % % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    % %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % % save(filename,'out','-append');

    % i=1; j=4; opt.alphaStep='FISTA_ADMM_NNL1';
    % opt.adaptiveStep=100; opt.cumuTol=10;
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');

    opt.skipIe=true;
    i=2; j=1; opt.alphaStep='FISTA_ADMM_NNL1';
    out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out','-append');

    return;

    % i=2; j=2; opt.alphaStep='NCG_PR';
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');

    % i=2; j=3; opt.alphaStep='IST_ADMM_NNL1';
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');
   
    opt.thresh=1e-6;
    i=2; j=4; opt.alphaStep='FISTA_ADMM_NNL1';
    out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out','-append');

    % opt.thresh=1e-6;
    % i=2; j=5; opt.alphaStep='FISTA_ADMM_NNL1';
    % opt.adaptiveStep=false;
    % out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
    %     conf.Psi,conf.Psit,conf.y,initSig,opt);
    % save(filename,'out','-append');
    % opt.adaptiveStep=true;

    opt.thresh=1e-6; opt.adaptiveStep=100; opt.cumuTol=10;
    i=2; j=6; opt.alphaStep='FISTA_ADMM_NNL1';
    out{i,j}=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out','-append');
end

if(any(runList==921))
    conf=ConfigCT('castSim','CircleMask','gpuPrj');
    conf.prjFull = 360; conf.prjNum = conf.prjFull;
    opt=conf.setup();
    y = conf.Phi(opt.trueAlpha); % equivalent to linear projection
    forSave=[conf.y y];
    polyy=linspace(min(conf.y),max(conf.y),1000);
    [y1,idx]=sort(conf.y); y2=y(idx);
    [y1,ia,ic]=unique(y1); y2=y2(ia);
    lineary=interp1(y1,y2,polyy,'spline');
    forsave=[polyy(:); lineary(:)];
    save('linearization.data','forSave','-ascii');

    return;

    load([mfilename '_021.mat']);
    fprintf('for non skiped Ie\n');
    t=1; noSkiped(:,t)=1:2000;
    for i=[1,3]
        t=t+1; noSkiped(1:length(out{1,i}.cost),t)=out{1,i}.cost;
        t=t+1; noSkiped(1:length(out{1,i}.RMSE),t)=out{1,i}.RMSE;
        t=t+1; noSkiped(1:length(out{1,i}.time),t)=out{1,i}.time;
    end
    mincost=reshape(noSkiped(:,[2,5]),[],1); mincost=min(mincost(mincost>0));
    noSkiped(:,[2,5])=noSkiped(:,[2,5])-mincost;
    save('costRmseTime.data','noSkiped','-ascii');
    return;
    figure; subplot(2,1,1);
    semilogy(out{1,1}.cost-mincost,'r'); hold on;
    semilogy(out{1,2}.cost-mincost,'g');
    semilogy(out{1,3}.cost-mincost,'b'); subplot(2,1,2);
    semilogy(out{1,1}.RMSE,'r'); hold on;
    semilogy(out{1,2}.RMSE,'g');
    semilogy(out{1,3}.RMSE,'b');
    figure; subplot(2,1,1);
    semilogy(out{1,1}.time, out{1,1}.cost-mincost,'r'); hold on;
    semilogy(out{1,2}.time, out{1,2}.cost-mincost,'g');
    semilogy(out{1,3}.time, out{1,3}.cost-mincost,'b'); subplot(2,1,2);
    semilogy(out{1,1}.time, out{1,1}.RMSE,'r'); hold on;
    semilogy(out{1,2}.time, out{1,2}.RMSE,'g');
    semilogy(out{1,3}.time, out{1,3}.RMSE,'b');

    return;
    t=1; skipped(:,t)=1:2000;
    for i=1:3
        out=out021{2,i};
        t=t+1; skipped(1:length(out.cost),t)=out.cost;
        t=t+1; skipped(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; skipped(1:length(out.time),t)=out.time;
    end
    mincost=reshape(skipped(:,[2,5,8]),[],1); mincost=min(mincost(mincost>0));
    skipped(:,[2,5,8])=skipped(:,[2,5,8])-mincost;
    save('costRmseTime_fixedI.data','skipped','-ascii');
end

% compare the effect forcing center model of spectrum 
if(any(runList==022))
    filename = [mfilename '_022.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    conf=ConfigCT();
    conf.PhiMode='gpuPrj';
    conf.prjFull = 360; conf.prjNum = conf.prjFull/2; opt.u = 1e-4;
    opt.debugLevel=6;
    opt=conf.setup(opt); initSig = maskFunc(conf.FBP(conf.y),opt.mask~=0);
    opt.maxIeSteps=1; opt.thresh=1e-6;

    opt.alphaStep='FISTA_ADMM_NNL1'; opt.CenterB=false; 
    out1=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out1','-append');

    opt.alphaStep='FISTA_ADMM_NNL1'; opt.CenterB=true;
    opt.correctCenterB=false;
    out2=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out2','-append');

    opt.alphaStep='FISTA_ADMM_NNL1'; opt.CenterB=true;
    opt.correctCenterB=true;
    out3=beamhardenSpline(conf.Phi,conf.Phit,...
        conf.Psi,conf.Psit,conf.y,initSig,opt);
    save(filename,'out3','-append');
end

if(any(runList==922))
    load([mfilename '_022.mat']);
    polymodel = Spline(out1.opt.spectBasis,out1.kappa);
    polymodel.setPlot(out1.opt.trueKappa,out1.opt.trueIota,out1.opt.epsilon);
    [tK,tU,kappa,Ie]=polymodel.plotSpectrum(out1.Ie);
    idx=find(tU==max(tU));
    q=tK(idx)/1; q=1/1.15
    tK=tK/q; tU=tU*q;
    i=1; forSave=[tK, tU];

    i=i+2; forSave(1:length(out1.kappa),i:i+1)=[out1.kappa,out1.Ie];
    figure; semilogx(tK,tU,'r'); hold on; loglog(out1.kappa,out1.Ie,'g-*');

    i=i+2; forSave(1:length(out2.kappa),i:i+1)=[out2.kappa,out2.Ie];
    loglog(out2.kappa,out2.Ie,'b-<');

    i=i+2; forSave(1:length(out3.kappa),i:i+1)=[out3.kappa,out3.Ie];
    loglog(out3.kappa,out3.Ie,'c-s');
    save('effectiveCenterB.data','forSave','-ascii');
    system(['mv effectiveCenterB.data ' paperDir]);
end

if(any(runList==122)) % dis, max AS step,
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

