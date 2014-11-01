function out = lasso(Phi,Phit,Psi,Psit,y,xInit,opt)
%lasso    Solve a linear problem
%
%                          0.5*||Φα-y||^2 + u*||Ψ'α||_1
%
%   where x can be anything, if opt.nu is set to be positive, x is assumed to be nonnegative.
%
%   Phi         The projection matrix implementation function handle
%   Phit        Transpose of Phi
%   Psi         Inverse wavelet transform matrix from wavelet coefficients
%               to image.
%   Psit        Transpose of Psi
%   y           Log scale of CT measurement y=-log(I^{mea}/I_0),
%   xInit       Initial value for the algorithm
%   opt         Structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%

if(~isfield(opt,'alphaStep')) opt.alphaStep='FISTA_L1'; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.8; end
if(~isfield(opt,'initStep')) opt.initStep='BB'; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=100; end
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
            computError = @(xxx) sqrNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
        case 2
            trueAlphaNorm=pNorm(opt.trueAlpha);
            computError = @(xxx) pNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
    end
end

if(opt.debugLevel>=2) figCost=1000; figure(figCost); end
if(opt.debugLevel>=3) figRes=1001; figure(figRes); end
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
        alphaStep=SpaRSA(2,alpha,1,opt.stepShrnk,Psi,Psit,alphaStep.M);
    case lower('FISTA_L1')
        alphaStep=FISTA_L1(1,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('FISTA_NN')
        alphaStep=FISTA_NN(2,alpha,1,opt.stepShrnk);
    case lower('FISTA_NNL1')
        alphaStep=FISTA_NNL1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('FISTA_ADMM_NNL1')
        alphaStep=FISTA_ADMM_NNL1(1,alpha,1,opt.stepShrnk,Psi,Psit);
        if(isfield(opt,'admmAbsTol')) alphaStep.admmAbsTol=opt.admmAbsTol; end
        if(isfield(opt,'admmTol')) alphaStep.admmTol=opt.admmTol; end
    case lower('IST_ADMM_NNL1')
        alphaStep=IST_ADMM_NNL1(1,alpha,1,opt.stepShrnk,Psi,Psit);
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
    case 'poisson'
        alphaStep.fArray{1} = @(aaa) Utils.poissonModel(aaa,Phi,Phit,y);
    case 'gaussian'
        alphaStep.fArray{1} = @(aaa) Utils.linearModel(aaa,Phi,Phit,y);
end
alphaStep.fArray{2} = @Utils.nonnegPen;
alphaStep.coef(1:2) = [1; opt.nu;];
if(any(strcmp(properties(alphaStep),'restart')) && (~opt.restart))
    alphaStep.restart=-1;
end
if(any(strcmp(properties(alphaStep),'adaptiveStep'))...
        && isfield(opt,'adaptiveStep'))
    alphaStep.adaptiveStep=opt.adaptiveStep;
end

if(opt.continuation || opt.fullcont)
    contIdx=1;
    if(length(opt.u)>1)
        alphaStep.u = opt.u(contIdx);
    else
        [~,g]=alphaStep.fArray{1}(alpha);
        alphaStep.u = opt.contEta*pNorm(Psit(g),inf);
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
    
    alphaStep.main();

    % if(p>=641) keyboard; end

    out.fVal(p,:) = (alphaStep.fVal(:))';
    out.cost(p) = alphaStep.cost;
    out.alphaSearch(p) = alphaStep.ppp;
    out.stepSize(p) = alphaStep.stepSize;
    if(opt.restart) out.restart(p)=alphaStep.restart; end
    if(collectNonInc) out.nonInc(p)=alphaStep.nonInc; end
    if(collectInnerSearch) out.innerSearch(p)=alphaStep.innerSearch; end;
    if(isfield(opt,'getBB') && opt.getBB)
        out.BB(p,:)=alphaStep.stepSizeInit('BB');
    end
    
    out.difAlpha(p)=relativeDif(alphaStep.alpha,alpha);
    if(p>1) out.difCost(p)=abs(out.cost(p)-out.cost(p-1))/out.cost(p); end

    alpha = alphaStep.alpha;

    str=sprintf([str ' cost=%-6g'],out.cost(p));

    if(isfield(opt,'trueAlpha'))
        out.RMSE(p)=computError(alpha);
        str=sprintf([str ' RSE=%g'],out.RMSE(p));
    end


    if(opt.continuation || opt.fullcont)
        out.uRecord(p,:)=[opt.u(end),alphaStep.u];
        str=sprintf([str ' u=%-6g'],alphaStep.u);
        temp=alphaStep.u/opt.u(end);
        if(opt.continuation)
            temp1=(opt.thresh*qThresh^(log(temp)/lnQU));
            temp1=max(temp1,1e-5);
        else
            temp1=opt.thresh;
        end
        out.contThresh(p)=temp1;
        if(temp>1 && out.difAlpha(p) < temp1 )
            out.contAlpha(:,contIdx)=alpha;
            if(isfield(opt,'trueAlpha')) out.contRMSE(contIdx)=out.RMSE(p); end
            contIdx=contIdx+1;
            if(length(opt.u)>1)
                alphaStep.u = opt.u(contIdx);
            else
                alphaStep.u = min(alphaStep.u*opt.contShrnk,opt.contEta*pNorm(Psit(alphaStep.grad),inf));
                alphaStep.u = max(alphaStep.u,opt.u);
            end
            alphaStep.reset();
        end
    end

    str=sprintf([str ' difAlpha=%g aSearch=%d'],out.difAlpha(p),alphaStep.ppp);
    if(p>1)
        str=sprintf([str ' difCost=%g'], out.difCost(p));
    end
    
    if(p>1 && opt.debugLevel>=2)
        set(0,'CurrentFigure',figCost);
        if(isfield(opt,'trueAlpha')) subplot(2,1,1); end
        semilogy(p-1:p,out.cost(p-1:p),'k'); hold on;
        title(sprintf('cost(%d)=%g',p,out.cost(p)));

        if(isfield(opt,'trueAlpha'))
            subplot(2,1,2);
            semilogy(p-1:p,out.RMSE(p-1:p)); hold on;
            title(sprintf('RMSE(%d)=%g',p,out.RMSE(p)));
        end
        drawnow;
    end

    if(length(out.fVal(p,:))>=1 && p>1 && opt.debugLevel>=3)
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
fprintf('\n');

end

