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
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'showImg')) opt.showImg=false; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
if(~isfield(opt,'continuation')) opt.continuation=false; end
if(~isfield(opt,'contShrnk')) opt.contShrnk=0.98; end
if(~isfield(opt,'contCrtrn')) opt.contCrtrn=1e-4; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-12; end
if(~isfield(opt,'maxItr')) opt.maxItr=1.5e3; end
% default to not use nonnegative constraints.
if(~isfield(opt,'nu')) opt.nu=0; end
if(~isfield(opt,'u')) opt.u=1e-4; end
if(~isfield(opt,'uMode')) opt.uMode='abs'; end
if(~isfield(opt,'muLustig')) opt.muLustig=1e-12; end
if(~isfield(opt,'errorType')) opt.errorType=1; end
if(~isfield(opt,'restart')) opt.restart=true; end
if(~isfield(opt,'noiseType')) opt.noiseType='gaussian'; end

alpha=xInit(:);

if(isfield(opt,'trueAlpha'))
    trueAlphaNorm=norm(opt.trueAlpha);
    switch opt.errorType
        case 0
            trueAlpha = opt.trueAlpha/trueAlphaNorm;
            computError= @(xxx) 1-(xxx(:)'*trueAlpha/norm(xxx(:)))^2;
        case 1
            computError = @(xxx) (norm(xxx(:)-opt.trueAlpha)/trueAlphaNorm)^2;
    end
end

if(opt.showImg && opt.debugLevel>=2) figCost=1000; figure(figCost); end
if(opt.showImg && opt.debugLevel>=3) figRes=1001; figure(figRes); end
if(opt.showImg && opt.debugLevel>=6) figAlpha=1002; figure(figAlpha); end

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
        alphaStep=FISTA_L1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('FISTA_NN')
        alphaStep=FISTA_NN(2,alpha,1,opt.stepShrnk);
    case lower('FISTA_NNL1')
        alphaStep=FISTA_NNL1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('FISTA_ADMM_NNL1')
        alphaStep=FISTA_ADMM_NNL1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case lower('IST_ADMM_NNL1')
        alphaStep=IST_ADMM_NNL1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_NNL1')}
        alphaStep=ADMM_NNL1(1,alpha,1,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_L1')}
        alphaStep = ADMM_L1(2,alpha,1,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_NN')}
        alphaStep = ADMM_NN(2,alpha,1,opt.stepShrnk,Psi,Psit);
end
if(strcmpi(opt.noiseType,'gaussian'))
    alphaStep.fArray{1} = @(aaa) Utils.linearModel(aaa,Phi,Phit,y);
else
    alphaStep.fArray{1} = @(aaa) Utils.poissonModel(aaa,Phi,Phit,y);
end
alphaStep.fArray{2} = @Utils.nonnegPen;
alphaStep.coef(1:2) = [1; opt.nu;];
if(any(strcmp(properties(alphaStep),'restart')) && (~opt.restart))
    alphaStep.restart=-1;
end

if(strcmpi(opt.uMode,'relative'))
    opt.u=opt.a*max(abs(Psit(Phit(y))));
end
if(opt.continuation)
    alphaStep.u = 0.1*max(abs(Psit(Phit(y))));
    alphaStep.u = min(alphaStep.u,opt.u*1000);
else alphaStep.u = opt.u;
end

tic; p=0; str=''; strlen=0; convThresh=0;
%figure(123); figure(386);
while(true)
    p=p+1;
    str=sprintf(['p=%-4d'],p);
    
    alphaStep.main();

    out.fVal(p,:) = (alphaStep.fVal(:))';
    out.cost(p) = alphaStep.cost;
    out.alphaSearch(p) = alphaStep.ppp;

    out.difAlpha(p)=norm(alphaStep.alpha(:)-alpha(:))/norm(alpha);
    if(p>1) out.difCost(p)=abs(out.cost(p)-out.cost(p-1))/out.cost(p); end

    alpha = alphaStep.alpha;

    if(opt.continuation)
        out.uRecord(p,:)=[opt.u,alphaStep.u];
        str=sprintf([str ' u=%-6g'],alphaStep.u);
        alphaStep.u = max(alphaStep.u*opt.contShrnk,opt.u);
    end
    if(isfield(opt,'trueAlpha'))
        out.RMSE(p)=computError(alpha);
    end

    str=sprintf([str ' cost=%-6g RSE=%g',...
        ' difAlpha=%g aSearch=%d'],...
        out.cost(p),out.RMSE(p), out.difAlpha(p), ...
        alphaStep.ppp);
    if(p>1)
        str=sprintf([str ' difCost=%g'], out.difCost(p));
    end
    
    if(opt.showImg && p>1 && opt.debugLevel>=2)
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
    if(p>1 && out.difAlpha(p)<opt.thresh && (alphaStep.u==opt.u))
        convThresh=convThresh+1;
    end
    if(p >= opt.maxItr || convThresh>10)
        if(opt.debugLevel==0)
            fprintf('%s',str);
        end
        break;
    end
end
out.alpha=alpha; out.p=p; out.opt = opt;
out.grad=alphaStep.grad;
fprintf('\n');

end

