function out = solver(Phi,Phit,Psi,Psit,y,xInit,opt)
%solver    Solve a sparse regularized problem
%
%                           L(x) + u*||Ψ'x||_1                        (1)
%
%   where L(x) is the likelihood function based on the measurements y. We
%   provide a few options through opt.noiseType to configurate the
%   measurement model. One popular case is when opt.noiseType='gaussian',
%
%                      0.5*||Φx-y||^2 + u*||Ψ'x||_1                   (2)
%
%   where u is set through opt.u. The methods provided through option
%   opt.alphaStep includes:
%
%   NPG         Solves the above problem with additional constrainsts: x>=0
%   NPGs        Solves exactly the above problem without nonnegativity
%               constraints
%   PG          The same with NPG, but without Nesterov's acceleration, not
%               recommended to use.
%
%   Parameters
%   ==========
%   Phi(Φ)      The projection matrix implementation function handle
%   Phit        Transpose of Phi
%   Psi(Ψ)      Inverse wavelet transform matrix from wavelet coefficients
%               to image.
%   Psit        Transpose of Psi, need to have ΨΨ'=I
%   y           The measurements according to different models:
%               opt.noiseType='gaussian'
%                   Ey = Φx, with gaussian noise
%               opt.noiseType='poisson'
%                   Ey = Φx+b, with poisson noise, where b is known provided by opt.bb
%               opt.noiseType='poissonLogLink'
%                   Ey=I_0 exp(-Φx) with poisson noise, where I_0 is unknown
%               opt.noiseType='poissonLogLink0'
%                   Ey=I_0 exp(-Φx) with poisson noise, where I_0 is known
%   xInit       Initial value for estimation of x
%   opt         Structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%

if(~isfield(opt,'alphaStep')) opt.alphaStep='NPGs'; end
if(~isfield(opt,'preSteps')) opt.preSteps=2; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.8; end
if(~isfield(opt,'initStep')) opt.initStep='hessian'; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=100; end   % currently not used
% default to not use nonnegative constraints.
if(~isfield(opt,'nu')) opt.nu=0; end
if(~isfield(opt,'u')) opt.u=1e-4; end
if(~isfield(opt,'uMode')) opt.uMode='abs'; end
if(~isfield(opt,'muLustig')) opt.muLustig=1e-12; end
if(~isfield(opt,'errorType')) opt.errorType=1; end
if(~isfield(opt,'restart')) opt.restart=true; end
if(~isfield(opt,'noiseType')) opt.noiseType='gaussian'; end
% continuation setup
if(~isfield(opt,'continuation')) opt.continuation=false; end
if(~isfield(opt,'contShrnk')) opt.contShrnk=0.5; end
if(~isfield(opt,'contCrtrn')) opt.contCrtrn=1e-3; end
if(~isfield(opt,'contEta')) opt.contEta=1e-2; end
if(~isfield(opt,'contGamma')) opt.contGamma=1e4; end
% find rse vs a, this option if true will disable "continuation"
if(~isfield(opt,'fullcont')) opt.fullcont=false; end

if(opt.fullcont)
    opt.continuation=false;
end

alpha=xInit;

if(isfield(opt,'trueAlpha'))
    switch opt.errorType
        case 0
            trueAlpha = opt.trueAlpha/pNorm(opt.trueAlpha);
            computError= @(xxx) 1-(innerProd(xxx,trueAlpha)^2)/sqrNorm(xxx);
        case 1
            trueAlphaNorm=sqrNorm(opt.trueAlpha);
            if(trueAlphaNorm==0) trueAlphaNorm=eps; end
            computError = @(xxx) sqrNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
        case 2
            trueAlphaNorm=pNorm(opt.trueAlpha);
            if(trueAlphaNorm==0) trueAlphaNorm=eps; end
            computError = @(xxx) pNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
    end
end

if(opt.debugLevel>=3) figCost=1000; figure(figCost); end
if(opt.debugLevel>=4) figRes=1001; figure(figRes); end
if(opt.debugLevel>=6) figAlpha=1002; figure(figAlpha); end

switch lower(opt.alphaStep)
    case lower('NCG_PR')
        alphaStep = NCG_PR(3,alpha);
        if(isfield(opt,'muHuber'))
            fprintf('use huber approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) huber(aaa,opt.muHuber,Psi,Psit);
        end
        if(isfield(opt,'muLustig'))
            fprintf('use lustig approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) lustigL1(aaa,opt.muLustig,Psi,Psit);
        end
    case {lower('SpaRSA')}
        alphaStep=SpaRSA(2,alpha,1,opt.stepShrnk,Psi,Psit,opt.M);
    case lower('NPGs')
        alphaStep=NPGs(1,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('FISTA_NN')
        alphaStep=FISTA_NN(2,alpha,1,opt.stepShrnk);
    case lower('FISTA_NNL1')
        alphaStep=FISTA_NNL1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('NPG')
        alphaStep=NPG(1,alpha,1,opt.stepShrnk,Psi,Psit);
        if(strcmpi(opt.noiseType,'poisson'))
            alphaStep.forcePositive=true;
            % opt.alphaStep='PG';
            % alphaStep=PG(1,alpha,1,opt.stepShrnk,Psi,Psit);
        end
    case lower('PG')
        alphaStep=PG(1,alpha,1,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_NNL1')}
        alphaStep=ADMM_NNL1(1,alpha,1,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_L1')}
        alphaStep = ADMM_L1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_NN')}
        alphaStep = ADMM_NN(2,alpha,1,opt.stepShrnk,Psi,Psit);
end
switch lower(opt.noiseType)
    case lower('poissonLogLink')
        alphaStep.fArray{1} = @(aaa) Utils.poissonModelLogLink(aaa,Phi,Phit,y);
    case lower('poissonLogLink0')
        alphaStep.fArray{1} = @(aaa) Utils.poissonModelLogLink0(aaa,Phi,Phit,y,opt.I0);
    case {'poisson', 'poisson2'} % poisson2 is for test
        if(isfield(opt,'bb'))
            temp=reshape(opt.bb,size(y));
        else
            temp=0;
        end
        alphaStep.fArray{1} = @(aaa) Utils.poissonModelAppr(aaa,Phi,Phit,y,temp);
        trueCost = @(aaa) Utils.poissonModelAppr(aaa,Phi,Phit,y,temp,1e-100);
    case 'gaussian'
        alphaStep.fArray{1} = @(aaa) Utils.linearModel(aaa,Phi,Phit,y);
end
alphaStep.fArray{2} = @Utils.nonnegPen;
alphaStep.coef(1:2) = [1; opt.nu;];
if(any(strcmp(properties(alphaStep),'restart')))
    if(~opt.restart) alphaStep.restart=-1; end
else
    opt.restart=false;
end
if(any(strcmp(properties(alphaStep),'adaptiveStep'))...
        && isfield(opt,'adaptiveStep'))
    alphaStep.adaptiveStep=opt.adaptiveStep;
end
if(any(strcmp(properties(alphaStep),'admmAbsTol'))...
        && isfield(opt,'admmAbsTol'))
    alphaStep.admmAbsTol=opt.admmAbsTol;
end
if(any(strcmp(properties(alphaStep),'admmTol'))...
        && isfield(opt,'admmTol'))
    alphaStep.admmTol=opt.admmTol;
end

if(opt.continuation || opt.fullcont)
    contIdx=1;
    inf_psit_grad=opt.u(1);
    if(length(opt.u)>1)
        alphaStep.u = opt.u(contIdx);
    else
        [~,g]=alphaStep.fArray{1}(alpha);
        inf_psit_grad=pNorm(Psit(g),inf);
        alphaStep.u = opt.contEta*inf_psit_grad;
        alphaStep.u = min(alphaStep.u,opt.u*opt.contGamma);
        alphaStep.u = max(alphaStep.u,opt.u);
        if(alphaStep.u*opt.contShrnk<=opt.u)
            opt.continuation=false;
            alphaStep.u=opt.u;
        end
        clear('g');
    end
    if(opt.continuation)
        qThresh = opt.contCrtrn/opt.thresh;
        lnQU = log(alphaStep.u/opt.u(end));
    end
else alphaStep.u = opt.u;
    fprintf('opt.u=%g\n',opt.u);
end

if(strcmpi(opt.initStep,'fixed'))
    alphaStep.stepSizeInit(opt.initStep,opt.L);
else alphaStep.stepSizeInit(opt.initStep);
end

if(any(strcmp(properties(alphaStep),'cumuTol'))...
        && isfield(opt,'cumuTol'))
    alphaStep.cumuTol=opt.cumuTol;
end

if(any(strcmp(properties(alphaStep),'innerSearch')))
    collectInnerSearch=true;
else
    collectInnerSearch=false;
end

if(any(strcmp(properties(alphaStep),'debug')))
    collectDebug=true;
    out.debug={};
else
    collectDebug=false;
end

if(any(strcmp(properties(alphaStep),'preSteps')))
    alphaStep.preSteps=opt.preSteps;
end

if(any(strcmp(properties(alphaStep),'nonInc')))
    collectNonInc=true;
else
    collectNonInc=false;
end

tic; p=0; strlen=0; convThresh=0;
%figure(123); figure(386);
while(true)
    p=p+1;
    str=sprintf('p=%-4d',p);
    if(p<=opt.preSteps && ~strcmpi(opt.initStep,'fixed'))
        temp=alphaStep.stepSizeInit(opt.initStep);
        alphaStep.t=min(alphaStep.t,temp);
    end
    
    %if(p>273) keyboard; end
    alphaStep.main();

    out.fVal(p,:) = (alphaStep.fVal(:))';
    out.cost(p) = alphaStep.cost;

    out.alphaSearch(p) = alphaStep.ppp;
    out.stepSize(p) = alphaStep.stepSize;
    if(opt.restart) out.restart(p)=alphaStep.restart; end
    if(collectNonInc) out.nonInc(p)=alphaStep.nonInc; end
    if(collectInnerSearch) out.innerSearch(p)=alphaStep.innerSearch; end;
    if(collectDebug && ~isempty(alphaStep.debug))
        out.debug{size(out.debug,1)+1,1}=p;
        out.debug{size(out.debug,1),2}=alphaStep.debug;
    end;
    if(opt.debugLevel>1)
        out.BB(p,1)=alphaStep.stepSizeInit('BB');
        out.BB(p,2)=alphaStep.stepSizeInit('hessian');
        % alphaStep.stepSizeInit('hessian',alpha);
    end
    
    out.difAlpha(p)=relativeDif(alphaStep.alpha,alpha);
    if(p>1) out.difCost(p)=abs(out.cost(p)-out.cost(p-1))/out.cost(p); end

    alpha = alphaStep.alpha;
    if(isfield(opt,'saveTrueCost') && opt.saveTrueCost && strcmpi(opt.noiseType,'poisson') )
            out.trueCost(p)=trueCost(alpha)+alphaStep.u*pNorm(Psit(alpha),1);
    end

    str=sprintf([str ' cost=%-6g'],out.cost(p));

    if(isfield(opt,'trueAlpha'))
        out.RMSE(p)=computError(alpha);
        str=sprintf([str ' RSE=%g'],out.RMSE(p));
    end


    if(opt.continuation || opt.fullcont)
        out.uRecord(p,:)=[opt.u(end),alphaStep.u,inf_psit_grad];
        str=sprintf([str ' u=%-6g'],alphaStep.u);
        temp=alphaStep.u/opt.u(end);
        if(opt.continuation)
            temp1=(opt.thresh*qThresh^(log(temp)/lnQU));
            temp1=max(temp1,opt.thresh*10);
        else
            temp1=opt.thresh;
        end
        if(temp>1) out.contThresh(p)=temp1; else out.contThresh(p)=opt.thresh; end;
        if(temp>1 && out.difAlpha(p) < temp1 )
            out.contAlpha{contIdx}=alpha;
            if(isfield(opt,'trueAlpha')) out.contRMSE(contIdx)=out.RMSE(p); end
            inf_psit_grad=pNorm(Psit(alphaStep.grad),inf);
            contIdx=contIdx+1;
            if(length(opt.u)>1)
                alphaStep.u = opt.u(contIdx);
            else
                alphaStep.u = min(alphaStep.u*opt.contShrnk,inf_psit_grad*opt.contEta);
                alphaStep.u = max(alphaStep.u,opt.u);
            end
            alphaStep.reset();
        end
    end

    str=sprintf([str ' difAlpha=%g aSearch=%d'],out.difAlpha(p),alphaStep.ppp);
    if(p>1)
        str=sprintf([str ' difCost=%g'], out.difCost(p));
    end
    
    if(p>1 && opt.debugLevel>=3)
        set(0,'CurrentFigure',figCost);
        if(isfield(opt,'trueAlpha')) subplot(2,1,1); end
        if(out.cost(p)>0)
            semilogy(p-1:p,out.cost(p-1:p),'k'); hold on;
            title(sprintf('cost(%d)=%g',p,out.cost(p)));
        end

        if(isfield(opt,'trueAlpha'))
            subplot(2,1,2);
            semilogy(p-1:p,out.RMSE(p-1:p)); hold on;
            title(sprintf('RMSE(%d)=%g',p,out.RMSE(p)));
        end
        drawnow;
    end

    if(length(out.fVal(p,:))>=1 && p>1 && opt.debugLevel>=4)
        set(0,'CurrentFigure',figRes);
        style={'r','g','b'};
        for i=1:length(out.fVal(p,:))
            subplot(3,1,i);
            semilogy(p-1:p,out.fVal(p-1:p,i),style{i}); hold on;
        end
        drawnow;
    end

    if(opt.debugLevel>=6)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        drawnow;
    end
    %if(mod(p,100)==1 && p>100) save('snapshotFST.mat'); end
    if(opt.debugLevel>=1)
        if(alphaStep.warned)
            fprintf('%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
        if(mod(p,opt.verbose)==0)
            strlen=0; fprintf('\n');
        else strlen = length(str);
        end
    end
    out.time(p)=toc;
    if(p>1 && out.difAlpha(p)<=opt.thresh && (alphaStep.u==opt.u(end)))
        convThresh=convThresh+1;
    end
    if(p >= opt.maxItr || convThresh>2)
        if(opt.debugLevel==0) fprintf('%s',str); end
        break;
    end
end
out.alpha=alpha; out.p=p; out.opt = opt;
out.grad=alphaStep.grad;
out.date=datestr(now);
fprintf('\nTime used: %d, cost=%g',out.time(end),out.cost(end));
if(isfield(opt,'trueAlpha'))
    if(opt.fullcont)
        idx = min(find(out.contRMSE==min(out.contRMSE)));
        if(out.contRMSE(idx)<out.RMSE(end))
            fprintf(', u=%g, RMSE=%g\n',opt.u(idx),out.contRMSE(idx));
        else
            fprintf(', RMSE=%g\n',out.RMSE(end));
        end
    else
        fprintf(', RMSE=%g\n',out.RMSE(end));
    end
else
    fprintf('\n');
end

end

