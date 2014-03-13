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
%   $Revision: 0.1 $ $Date: Thu 13 Mar 2014 04:01:48 PM CDT
%

if(~isfield(opt,'alphaStep')) opt.alphaStep='FISTA_L1'; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'showImg')) opt.showImg=false; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=0; end
if(~isfield(opt,'continuation')) opt.continuation=false; end
if(~isfield(opt,'contShrnk')) opt.contShrnk=0.5; end
if(~isfield(opt,'contCrtrn')) opt.contCrtrn=0.5; end
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.thresh=1e3; end
% default to not use nonnegative constraints.
if(~isfield(opt,'nu')) opt.nu=0; end
if(~isfield(opt,'u')) opt.u=1e-4; end
if(~isfield(opt,'muLustig')) opt.muLustig=1e-12; end

alpha=xInit(:);

if(isfield(opt,'trueAlpha'))
    trueAlpha = opt.trueAlpha/norm(opt.trueAlpha);
end

if(opt.showImg)
    figRes=1000; figure(figRes);
    figAlpha=1001; figure(figAlpha);
else figRes=0; figAlpha=0;
end

%max(Imea./(exp(-atten(Phi,alpha)*mu')*Ie))
if(interiorPointAlpha)
    nonneg=@nonnegLogBarrier; 
else nonneg=@nonnegPen; 
end
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
        alphaStep=SpaRSA(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit,alphaStep.M);
    case lower('FISTA_L1')
        alphaStep=FISTA_L1(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        alphaStep.coef(1:2) = [1; 0];
        alphaStep.fArray{2} = nonneg;
        out.fVal=zeros(opt.maxItr,3);
    case {lower('ADMM_NNL1')}
        alphaStep=ADMM_NNL1(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        out.fVal=zeros(opt.maxItr,2);
    case {lower('ADMM_L1')}
        alphaStep = ADMM_L1(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        alphaStep.coef(1:2) = [1; 0];
        alphaStep.fArray{2} = nonneg;
        out.fVal=zeros(opt.maxItr,3);
end
alphaStep.fArray{1} = @(aaa) linearModel(aaa,Phi,Phit,y);
alphaStep.fArray{2} = nonneg;
alphaStep.stepShrnk = opt.stepShrnk;
alphaStep.coef(1:3) = [1; opt.nu; opt.u];

if(opt.continuation)
    alphaStep.u = 0.1*max(abs(Psit(Phit(y))));
else alphaStep.u = opt.u;
end

tic; p=0; str=''; strlen=0;
while(true)
    p=p+1;
    str=sprintf([str '\np=%d '],p);
    
    alphaStep.main();

    out.fVal(p,:) = (alphaStep.fVal(:))';
    out.cost(p) = alphaStep.cost;
    out.difAlpha(p) = norm(alphaStep.alpha(:)-alpha(:))^2;
    out.alphaSearch(p) = alphaStep.ppp;
    alpha = alphaStep.alpha;

    if(opt.continuation && p>1 && ...
            abs(out.cost(p)-out.cost(p-1))/out.cost(p)<opt.contCrtrn && ...
            alphaStep.u*opt.contShrnk>opt.u)
        alphaStep.u = max(alphaStep.u*opt.contShrnk,opt.u);
        fprintf('\nnew u= %g\n',alphaStep.u);
        alphaStep.warned = true;
    end
    %if(out.stepSz~=s1) fprintf('lineSearch is useful!!\n'); end
    if(isfield(opt,'trueAlpha'))
        out.RMSE(p)=1-(alpha'*trueAlpha/norm(alpha))^2;
    end

    str=sprintf([str 'cost=%-10g RSE=%-10g ',...
        'dAlpha=%-10g aSearch=%d '],...
        out.cost(p),out.RMSE(p), out.difAlpha(p), ...
        alphaStep.ppp);
    if(p>1)
        str=sprintf([str 'pdObjAlpha=%g%% '],...
            (out.cost(p-1)-out.cost(p))/out.cost(p)*100);
    end
    
    if(opt.showImg && p>1)
        set(0,'CurrentFigure',figRes);
        subplot(2,2,1);
        semilogy(p-1:p,out.fVal(p-1:p,1),'r'); hold on;
        semilogy(p-1:p,out.cost(p-1:p),'k');
        title(sprintf('cost(%d)=%g',p,out.cost(p)));

        subplot(2,2,2);
        if(length(out.fVal(p,:))>1)
            semilogy(p-1:p,out.fVal(p-1:p,2),'g');
        end
        if(length(out.fVal(p,:))>2)
            semilogy(p-1:p,out.fVal(p-1:p,3),'b');
        end

        if(isfield(opt,'trueAlpha'))
            subplot(2,1,2);
            semilogy(p-1:p,out.RMSE(p-1:p)); hold on;
            title(sprintf('RMSE(%d)=%g',p,out.RMSE(p)));
        end
        drawnow;
    end
    
    if(figAlpha)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        %showImgMask(Qmask-Qmask1/2,opt.mask);
        %title(['size of Q=' num2str(length(Q))]);
        if(~isempty(IeStep.zmf))
            title(sprintf('zmf=(%g,%g)', IeStep.zmf(1), IeStep.zmf(2)))
        end
        drawnow;
    end
    %if(mod(p,100)==1 && p>100) save('snapshotFST.mat'); end
    if(opt.visible && p>1)
        if(alphaStep.warned || IeStep.warned)
            fprintf('%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
        strlen = length(str);
        str='';
    end
    out.time(p)=toc;
    if(p >= opt.maxItr) break; end
    if(p>1 && abs(out.cost(p)-out.cost(p-1))/out.cost(p)<opt.thresh)
        break
    end
end

out.alpha=alpha; out.p=p; out.opt = opt;
fprintf('\n');

end

