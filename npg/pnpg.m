function out = pnpg(NLL,proximal,xInit,opt)
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
%               constraints;
%               Note that, NPGs support weighted l1 norm by specifying
%               option opt.weight
%   PG          The same with NPG, but without Nesterov's acceleration, not
%               recommend to use.
%
%   Parameters
%   ==========
%   Phi(Φ)      The projection matrix or its implementation function handle
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
%               opt.noiseType='logistic'
%                   Ey=exp(Φx+b)./(1+exp(Φx+b)) with Bernoulli noise, where b is optional
%   xInit       Initial value for estimation of x
%   opt         Structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Reference:
%   Author: Renliang Gu (gurenliang@gmail.com)

% default to not use any constraints.
if(~isfield(opt,'prj_C')) opt.prj_C=@(x)x; end
if(~isfield(opt,'u')) opt.u=1e-4; end

if(~isfield(opt,'alphaStep')) opt.alphaStep='NPGs'; end
if(~isfield(opt,'stepIncre')) opt.stepIncre=0.9; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'initStep')) opt.initStep='hessian'; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'saveXtrace')) opt.saveXtrace=false; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'errorType')) opt.errorType=1; end
if(~isfield(opt,'restart')) opt.restart=true; end
if(~isfield(opt,'gamma')) gamma=2; else gamma=opt.gamma; end
if(~isfield(opt,'b')) b=0.25; else b=opt.b; end

% For debug purpose
if(~isfield(opt,'collectNonInc'))
    collectNonInc=false;
else
    if(opt.collectNonInc) collectNonInc=true; end
end
if(~isfield(opt,'collectTheta'))
    collectTheta=false;
else
    if(opt.collectTheta) collectTheta=true; end
end

if(isfield(opt,'trueAlpha'))
    switch opt.errorType
        case 0
            trueAlpha = opt.trueAlpha/pNorm(opt.trueAlpha);
            computError= @(xxx) 1-(realInnerProd(xxx,trueAlpha)^2)/sqrNorm(xxx);
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

if(debug.level>=4)
    figCost=1000; figure(figCost);
    if(debug.level>=5)
        figRes=1001; figure(figRes);
        if(debug.level>=6)
            figAlpha=1002; figure(figAlpha);
        end
    end
end

switch lower(opt.alphaStep)
    case lower('NCG_PR')
    case {lower('NPGs')}
    case lower('PG')
    case {lower('NPG')}
    case {lower('PNPG')}
        oneStep=@pnpgStep;
    case lower('ATs')
    case lower('AT')
    case lower('GFB')
    case lower('Condat')
    case lower('FISTA_NN')
    case lower('FISTA_NNL1')
    case {lower('ADMM_NNL1')}
    case {lower('ADMM_L1')}
    case {lower('ADMM_NN')}
    otherwise

end

% specify the domain of the NLL via a projection operator
alphaStep.prj_C=opt.prj_C;
alpha=double(xInit);
debug=Debug(opt.debugLevel);

% determin whether to have convex set C constraints

% determin the Method to use

% determin the data fedelity term

if(strcmpi(opt.initStep,'fixed'))
    stepSizeInit(opt.initStep,opt.Lip);
else stepSizeInit(opt.initStep);
end

if(proximal.iterative && nargout(proximal.op)>1)
    collectInnerSearch=true;
else
    collectInnerSearch=false;
end

if(any(strcmp(properties(alphaStep),'debug')) && alphaStep.hasLog)
    collectDebug=true;
    out.debug={};
else
    collectDebug=false;
end

% print start information
if(debug.level>=1)
    fprintf('%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('??? decide whether to have "projected" Nestrov''s Proximal Gradient Method (%s) %s_%s',opt.proximal,opt.alphaStep,opt.noiseType);
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','Itr');
    str=sprintf([str ' %12s'],'Objective');
    if(isfield(opt,'trueAlpha'))
        str=sprintf([str ' %12s'], 'Error');
    end
    if(opt.continuation || opt.fullcont)
        str=sprintf([str ' %12s'],'u');
    end
    str=sprintf([str ' %12s %4s'], '|d α|/|α|', 'αSrh');
    str=sprintf([str ' %12s'], '|d Obj/Obj|');
    str=sprintf([str '\t u=%g'],opt.u);
    fprintf('%s\n%s',str,repmat( '-', 1, 80 ) );
end

tic; p=0; convThresh=0; theta=1;
%figure(123); figure(386);
while(true)
    p=p+1;
    %if(mod(p,100)==1 && p>100) save('snapshotFST.mat'); end

    oneStep;

    out.cost(p) = alphaStep.cost;
    out.alphaSearch(p) = alphaStep.ppp;
    out.stepSize(p) = alphaStep.stepSize;
    if(opt.restart) out.restart(p)=alphaStep.restart; end
    if(collectNonInc) out.nonInc(p)=alphaStep.nonInc; end
    if(collectTheta) out.theta(p)=alphaStep.theta; end
    if(collectInnerSearch) out.innerSearch(p)=alphaStep.innerSearch; end;
    if(collectDebug && ~isempty(alphaStep.debug))
        out.debug{size(out.debug,1)+1,1}=p;
        out.debug{size(out.debug,1),2}=alphaStep.debug.log;
    end;
    if(debug.level>1)
        out.BB(p,1)=alphaStep.stepSizeInit('BB');
        out.BB(p,2)=alphaStep.stepSizeInit('hessian');
    end

    out.difAlpha(p)=relativeDif(alphaStep.alpha,alpha);
    if(p>1) out.difCost(p)=abs(out.cost(p)-out.cost(p-1))/out.cost(p); end

    alpha = alphaStep.alpha;

    if(mod(p,opt.verbose)==1) debug.println(1); else debug.clear(1); end
    debug.print(1,sprintf(' %5d',p));
    debug.print(1,sprintf(' %12g',out.cost(p)));

    if(isfield(opt,'trueAlpha'))
        out.RMSE(p)=computError(alpha);
        debug.print(1,sprintf(' %12g',out.RMSE(p)));
    end

    if(opt.saveXtrace) out.alphaTrace(:,p)=alpha; end

    debug.print(1,sprintf(' %12g %4d',out.difAlpha(p),alphaStep.ppp));
    if(p>1)
        debug.print(1,sprintf(' %12g', out.difCost(p)));
    else
        debug.print(1,sprintf(' %12s', ' '));
    end

    if(p>1 && debug.level>=4)
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

    if(length(out.fVal(p,:))>=1 && p>1 && debug.level>=5)
        set(0,'CurrentFigure',figRes);
        style={'r','g','b'};
        for i=1:length(out.fVal(p,:))
            subplot(3,1,i);
            semilogy(p-1:p,out.fVal(p-1:p,i),style{i}); hold on;
        end
        drawnow;
    end

    if(debug.level>=6)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        drawnow;
    end

    out.time(p)=toc;

    if(p>1 && out.difAlpha(p)<=opt.thresh && (alphaStep.u==opt.u(end)))
        convThresh=convThresh+1;
    end

    if(p >= opt.maxItr || (convThresh>2 && p>opt.minItr))
        break;
    end
end
out.alpha=alpha; out.p=p; out.opt = opt;
out.grad=alphaStep.grad;
out.date=datestr(now);
if(debug.level>=0)
    fprintf('\nCPU Time: %g, objective=%g',out.time(end),out.cost(end));
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

function pnpgStep
    debug.clearLog();

    ppp=0; incStep=false; goodMM=true;
    if(adaptiveStep && cumu>=cumuTol)
        % adaptively increase the step size
        t=t*stepIncre;
        cumu=0;
        incStep=true;
    end
    % start of line Search
    while(true)
        ppp = ppp+1;

        B=t/preT;
        %newTheta=(1+sqrt(1+4*B*theta^2))/2;
        if(p==1)
            newTheta=1;
        else
            newTheta=1/gamma+sqrt(b+B*theta^2);
        end
        xbar=alpha+(theta -1)/newTheta*(alpha-preAlpha);
        xbar=prj_C(xbar);

        [oldCost,grad] = NLL(xbar);

        newX=proximal.op(xbar-grad/t,u/t,admmTol*difAlpha,maxInnerItr,...
            pInit);

        newCost=NLL(newX);
        if(majorizationHolds(newX-xbar,newCost,oldCost,[],grad,t))
            break;
        else
            if(ppp<=20 && t>0)
                t=t/stepShrnk;
                % Penalize if there is a step size increase just now
                if(incStep)
                    if(incCumuTol)
                        cumuTol=cumuTol+4;
                    end
                    incStep=false;
                end
            else  % don't know what to do, mark on debug and break
                if(t<0)
                    error('\n PNPG is having a negative step size, do nothing and return!!');
                end
                goodMM=false;
                debug.appendLog('_FalseMM');
                break;
            end
        end
    end
    stepSize = 1/t;
    newObj = newCost+u*proximal.penalty(newX);
    objBar = oldCost+u*proximal.penalty(xbar);

    if((newObj-cost)>0)
        if(goodMM && pNorm(xbar-alpha,1)~=0 && restart>=0) % if has monmentum term, restart
            theta=1;
            debug.appendLog('_Restart');
            debug.printWithoutDel(1,'\t restart');
            continue;
        elseif((~goodMM) || (objBar<newObj))
            if(~goodMM)
                debug.appendLog('_Reset');
                reset();
            end
            if(proximal.iterative)
                if(proximal.steps<maxInnerItr && admmTol>1e-6)
                    admmTol=admmTol/10;
                    debug.printWithoutDel(1,...
                        sprintf('\n decrease admmTol to %g',admmTol));
                    continue;
                    %% IMPORTANT! if not requir too high accuracy
                    %% use 1e3 for maxInnerItr
                elseif(proximal.steps>=maxInnerItr && maxInnerItr<maxPossibleInnerItr)
                    maxInnerItr=maxInnerItr*10;
                    debug.printWithoutDel(1,...
                        sprintf('\n increase maxInnerItr to %g',maxInnerItr));
                    continue;
                end
            end

            % give up and force it to converge
            debug.appendLog('_ForceConverge');
            newObj=cost;  newX=alpha;
            innerSearch=0;
        end
    else
        if(proximal.iterative)
            proximal.setInit();
            innerSearch=proximal.steps;
        else
            innerSearch=0;
        end
    end
    theta = newTheta; preAlpha = alpha;
    cost = newObj;
    difAlpha = relativeDif(alpha,newX);
    alpha = newX;
    preT=t;

    if(ppp==1 && adaptiveStep)
        cumu=cumu+1;
    else
        cumu=0;
    end
end

    function t=stepSizeInit(select,delta)
        switch (lower(select))
            case 'bb'   % use BB method to guess the initial stepSize
                if(~exist('delta','var'))
                    delta=1e-5;
                end
                [~,grad1] = func(alpha);
                temp = delta*grad1/pNorm(grad1);
                temp = alpha-prj_C(alpha-temp);
                [~,grad2] = func(alpha-temp);
                t = abs(realInnerProd(grad1-grad2,temp))/sqrNorm(temp);
            case 'hessian'
                [~,grad1,hessian] = func(alpha);
                if(isempty(hessian))
                    if(~exist('delta','var')) delta=1e-5; end
                    temp = delta*grad1/pNorm(grad1);
                    temp = alpha-prj_C(alpha-temp);
                    [~,grad2] = func(alpha-temp);
                    t = abs(realInnerProd(grad1-grad2,temp))/sqrNorm(temp);
                else
                    t = hessian(grad1,2)/sqrNorm(grad1);
                end
            case 'fixed'
                t = delta;
            otherwise
                error('unkown selection for initial step');
        end
        if(isnan(t)) t=ones(size(t)); end
        if(min(t)==0) keyboard;  t=eps; end;
        if(nargout==0)  keyboard; t=min(t); end
    end

function test = majorizationHolds(x_minus_y,fx,fy,dfx,dfy,L)
    % This function tests whether
    %      f(x) ≤ f(y)+(x-y)'*∇f(y)+ 0.5*L*||x-y||^2
    % holds.

    if(~isempty(dfx) && abs(fx-fy)/max(max(fx,fy),1) < 1e-10)
        % In this case, use stronger condition to avoid numerical issue
        test=(realInnerProd(x_minus_y,dfx-dfy) <= 0.5*L*sqrNorm(x_minus_y));
    else
        test=(fx-fy<=realInnerProd(x_minus_y,dfy)+0.5*L*sqrNorm(x_minus_y));
    end
end

