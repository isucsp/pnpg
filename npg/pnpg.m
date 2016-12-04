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
if(~isfield(opt,'prj_C') || isempty(opt.prj_C))
    opt.prj_C=@(x) x;
end

if(~isfield(opt,'u')) opt.u=1e-4; end

if(~isfield(opt,'alphaStep')) opt.alphaStep='PNPG'; end
if(~isfield(opt,'stepIncre')) opt.stepIncre=0.9; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'initStep')) opt.initStep='hessian'; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'saveXtrace')) opt.saveXtrace=false; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'Lip')) opt.Lip=nan; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'errorType')) opt.errorType=1; end
if(~isfield(opt,'restart')) opt.restart=true; end
if(~isfield(opt,'gamma')) gamma=2; else gamma=opt.gamma; end
if(~isfield(opt,'b')) b=0.25; else b=opt.b; end

if(~isfield(opt,'relInnerThresh')) opt.relInnerThresh=1e-2; end
if(~isfield(opt,'cumuTol')) opt.cumuTol=4; end
if(~isfield(opt,'incCumuTol')) opt.incCumuTol=true; end
if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=true; end
if(~isfield(opt,'maxInnerItr')) opt.maxInnerItr=100; end
if(~isfield(opt,'maxPossibleInnerItr')) opt.maxPossibleInnerItr=1e3; end

% Output options
if(~isfield(opt,'outDetail')) opt.outDetail=false; end
if(~isfield(opt,'NLL_Pen') || ~opt.outDetail) opt.NLL_Pen=false; end

% For debug purpose
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

debug=Debug(opt.debugLevel);
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
        if(opt.adaptiveStep)
            oneStep=@pnpgAdaptiveStep;
        else
            oneStep=@pnpgNonAdaptiveStep;
        end
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

% determin whether to have convex set C constraints

% determin the Method to use

% determin the data fedelity term

if(proximal.iterative && abs(nargout(proximal.op))>1)
    collectInnerItr=true;
else
    collectInnerItr=false;
end

if(nargin(proximal.op)==1)
    proximalOp=proximal.op;
    proximal.op=@(a,u) proximalOp(a);
end

collectDebug=true;
out.debug={};

% print start information
if(debug.level>=1)
    fprintf('%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Projected Nestrov''s Proximal-Gradient Method (%s)',opt.alphaStep);
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','Itr');
    str=sprintf([str ' %12s'],'Objective');
    if(isfield(opt,'trueAlpha'))
        str=sprintf([str ' %12s'], 'Error');
    end
    str=sprintf([str ' %12s %4s'], '|dx|/|x|', 'αSrh');
    str=sprintf([str ' %12s'], '|d Obj/Obj|');
    str=sprintf([str '\t u=%g'],opt.u);
    fprintf('%s\n%s',str,repmat( '-', 1, 80 ) );
end

t=stepSizeInit(opt.initStep,opt.Lip);

tic; itr=0; convThresh=0; theta=1;
restart=0;   % make this value negative to disable restart

% The variables that will be use in oneStep
% ??? cost can be more accurate initially
cost=inf; preT=t; preX=0; x=xInit; difX=1; nLineSearch=0;

if(opt.outDetail)
    innerItr=[]; stepSize=0;
    theta;
    grad=[];
    if(opt.NLL_Pen)
        NLLVal=0; penVal=0;
    end
end


if(proximal.iterative)
    pInit=[];
end

if(opt.adaptiveStep)
    cumu=0;
end

while(true)
    itr=itr+1;
    %if(mod(itr,100)==1 && itr>100) save('snapshotFST.mat'); end

    oneStep();

    out.cost(itr) = cost;
    if(opt.outDetail)
        out.alphaSearch(itr) = nLineSearch;
        out.stepSize(itr) = stepSize;
        if(opt.NLL_Pen)
            out.NLLVal(itr)=NLLVal;
            out.penVal(itr)=penVal;
        end
    end
    if(opt.restart) out.restart(itr)=restart; end
    if(collectTheta) out.theta(itr)=theta; end
    if(collectInnerItr) out.innerItr(itr)=innerItr; end;
    if(collectDebug && ~isempty(debug.log))
        out.debug{size(out.debug,1)+1,1}=itr;
        out.debug{size(out.debug,1),2}=debug.log;
    end;
    if(debug.level>1)
        out.BB(itr,1)=stepSizeInit('BB');
        out.BB(itr,2)=stepSizeInit('hessian');
    end

    out.difX(itr)=difX;
    if(itr>1) out.difCost(itr)=abs(out.cost(itr)-out.cost(itr-1))/out.cost(itr); end

    if(mod(itr,opt.verbose)==1) debug.println(1); else debug.clear(1); end
    debug.print(1,sprintf(' %5d',itr));
    debug.print(1,sprintf(' %12g',out.cost(itr)));

    if(isfield(opt,'trueAlpha'))
        out.RMSE(itr)=computError(x);
        debug.print(1,sprintf(' %12g',out.RMSE(itr)));
    end

    if(opt.saveXtrace) out.alphaTrace(:,itr)=x; end

    debug.print(1,sprintf(' %12g %4d',out.difX(itr),nLineSearch));
    if(itr>1)
        debug.print(1,sprintf(' %12g', out.difCost(itr)));
    else
        debug.print(1,sprintf(' %12s', ' '));
    end

    if(itr>1 && debug.level>=4)
        set(0,'CurrentFigure',figCost);
        if(isfield(opt,'trueAlpha')) subplot(2,1,1); end
        if(out.cost(itr)>0)
            semilogy(itr-1:itr,out.cost(itr-1:itr),'k'); hold on;
            title(sprintf('cost(%d)=%g',itr,out.cost(itr)));
        end

        if(isfield(opt,'trueAlpha'))
            subplot(2,1,2);
            semilogy(itr-1:itr,out.RMSE(itr-1:itr)); hold on;
            title(sprintf('RMSE(%d)=%g',itr,out.RMSE(itr)));
        end
        drawnow;
    end

    if(opt.NLL_Pen && itr>1 && debug.level>=5)
        set(0,'CurrentFigure',figRes);
        subplot(2,1,1);
        semilogy(itr-1:itr,out.NLLVal(itr-1:itr),'r'); hold on;
        subplot(2,1,2);
        semilogy(itr-1:itr,out.penVal(itr-1:itr),'b'); hold on;
        drawnow;
    end

    if(debug.level>=6)
        set(0,'CurrentFigure',figAlpha); showImgMask(x,opt.mask);
        drawnow;
    end

    out.time(itr)=toc;

    if(itr>1 && out.difX(itr)<=opt.thresh )
        convThresh=convThresh+1;
    end

    if(itr >= opt.maxItr || (convThresh>2 && itr>opt.minItr))
        break;
    end
end
out.x=x; out.itr=itr; out.opt = opt;
if(opt.outDetail)
    out.grad=grad;
end
out.date=datestr(now);
if(debug.level>=0)
    fprintf('\nCPU Time: %g, objective=%g',out.time(end),out.cost(end));
    if(isfield(opt,'trueAlpha'))
        fprintf(', RMSE=%g\n',out.RMSE(end));
    else
        fprintf('\n');
    end
end

function pnpgAdaptiveStep
    debug.clearLog();

    while(true)
        nLineSearch=0; incStep=false; goodMM=true;
        if(cumu>=opt.cumuTol)
            % adaptively increase the step size
            t=t*opt.stepIncre;
            cumu=0;
            incStep=true;
        end
        % start of line Search
        while(true)
            nLineSearch = nLineSearch+1;

            B=t/preT;
            %newTheta=(1+sqrt(1+4*B*theta^2))/2;
            if(itr==1)
                newTheta=1;
            else
                newTheta=1/gamma+sqrt(b+B*theta^2);
            end
            xbar=x+(theta-1)/newTheta*(x-preX);
            xbar=opt.prj_C(xbar);

            [oldCost,grad] = NLL(xbar);

            if(proximal.iterative)
                [newX,innerItr_,pInit_]=proximal.op(xbar-grad/t,opt.u/t,opt.relInnerThresh*difX,opt.maxInnerItr,...
                    pInit);
            else
                newX=proximal.op(xbar-grad/t,opt.u/t);
            end

            newCost=NLL(newX);
            if(majorizationHolds(newX-xbar,newCost,oldCost,[],grad,t))
                break;
            else
                if(nLineSearch<=20 && t>0)
                    t=t/opt.stepShrnk;
                    % Penalize if there is a step size increase just now
                    if(incStep)
                        if(opt.incCumuTol)
                            opt.cumuTol=opt.cumuTol+4;
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
        newPen = proximal.penalty();
        newObj = newCost+opt.u*newPen;

        if((newObj-cost)>0)
            oldPen = proximal.penalty(xbar);
            objBar = oldCost+opt.u*oldPen;
            if(goodMM && pNorm(xbar-x,1)~=0 && restart>=0) % if has monmentum term, restart
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
                    if(proximal.steps<opt.maxInnerItr && opt.relInnerThresh>1e-6)
                        opt.relInnerThresh=opt.relInnerThresh/10;
                        debug.printWithoutDel(1,...
                            sprintf('\n decrease relInnerThresh to %g',...
                            opt.relInnerThresh));
                        continue;
                        %% IMPORTANT! if not requir too high accuracy
                        %% use 1e3 for maxInnerItr
                    elseif(proximal.steps>=opt.maxInnerItr &&...
                            opt.maxInnerItr<opt.maxPossibleInnerItr)
                        opt.maxInnerItr=opt.maxInnerItr*10;
                        debug.printWithoutDel(1,...
                            sprintf('\n increase maxInnerItr to %g',opt.maxInnerItr));
                        continue;
                    end
                end

                % give up and force it to converge
                debug.appendLog('_ForceConverge');
                newObj=cost;  newX=x;
                innerItr=0;
            end
        else
            if(proximal.iterative)
                pInit=pInit_;
                innerItr=innerItr_;
            else
                innerItr=0;
            end
        end
        theta = newTheta; preX = x;
        cost = newObj;
        difX = relativeDif(x,newX);
        x = newX;
        preT=t;
        penVal=newPen;
        NLLVal=newCost;

        if(nLineSearch==1)
            cumu=cumu+1;
        else
            cumu=0;
        end

        break;
    end
    function reset()
        theta=1; t=min([t;max(stepSizeInit('hessian'))]);
    end
end

function pnpgNonAdaptiveStep
    debug.clearLog();

    while(true)
        nLineSearch=0; goodMM=true;

        %newTheta=(1+sqrt(1+4*B*theta^2))/2;
        if(itr==1)
            newTheta=1;
        else
            newTheta=1/gamma+sqrt(b+theta^2);
        end
        xbar=x+(theta-1)/newTheta*(x-preX);
        xbar=opt.prj_C(xbar);
        [oldCost,grad] = NLL(xbar);

        % start of line Search
        while(true)
            nLineSearch = nLineSearch+1;

            if(proximal.iterative)
                [newX,innerItr_,pInit_]=proximal.op(xbar-grad/t,opt.u/t,opt.relInnerThresh*difX,opt.maxInnerItr,...
                    pInit);
            else
                newX=proximal.op(xbar-grad/t,opt.u/t);
            end

            newCost=NLL(newX);
            if(majorizationHolds(newX-xbar,newCost,oldCost,[],grad,t))
                break;
            else
                if(nLineSearch<=20 && t>0)
                    t=t/opt.stepShrnk;
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
        newPen = proximal.penalty();
        newObj = newCost+opt.u*newPen;

        if((newObj-cost)>0)
            oldPen = proximal.penalty(xbar);
            objBar = oldCost+opt.u*oldPen;
            if(goodMM && pNorm(xbar-x,1)~=0 && restart>=0) % if has monmentum term, restart
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
                    if(proximal.steps<opt.maxInnerItr && opt.relInnerThresh>1e-6)
                        opt.relInnerThresh=opt.relInnerThresh/10;
                        debug.printWithoutDel(1,...
                            sprintf('\n decrease relInnerThresh to %g',...
                            opt.relInnerThresh));
                        continue;
                        %% IMPORTANT! if not requir too high accuracy
                        %% use 1e3 for maxInnerItr
                    elseif(proximal.steps>=opt.maxInnerItr &&...
                            opt.maxInnerItr<opt.maxPossibleInnerItr)
                        opt.maxInnerItr=opt.maxInnerItr*10;
                        debug.printWithoutDel(1,...
                            sprintf('\n increase maxInnerItr to %g',opt.maxInnerItr));
                        continue;
                    end
                end

                % give up and force it to converge
                debug.appendLog('_ForceConverge');
                newObj=cost;  newX=x;
                innerItr=0;
            end
        else
            if(proximal.iterative)
                pInit=pInit_;
                innerItr=innerItr_;
            else
                innerItr=0;
            end
        end
        theta = newTheta; preX = x;
        cost = newObj;
        difX = relativeDif(x,newX);
        x = newX;
        penVal=newPen;
        NLLVal=newCost;
        break;
    end
    function reset()
        theta=1; t=min([t;max(stepSizeInit('hessian'))]);
    end
end

function t=stepSizeInit(select,Lip,delta)
    switch (lower(select))
        case 'bb'   % use BB method to guess the initial stepSize
            if(~exist('delta','var'))
                delta=1e-5;
            end
            [~,grad1] = NLL(x);
            temp = delta*grad1/pNorm(grad1);
            temp = x-opt.prj_C(x-temp);
            [~,grad2] = NLL(x-temp);
            t = abs(realInnerProd(grad1-grad2,temp))/sqrNorm(temp);
        case 'hessian'
            [~,grad1,hessian] = NLL(x);
            if(isempty(hessian))
                if(~exist('delta','var')) delta=1e-5; end
                temp = delta*grad1/pNorm(grad1);
                temp = x-opt.prj_C(x-temp);
                [~,grad2] = NLL(x-temp);
                t = abs(realInnerProd(grad1-grad2,temp))/sqrNorm(temp);
            else
                t = hessian(grad1,2)/sqrNorm(grad1);
            end
        case 'fixed'
            t = Lip;
        otherwise
            error('unkown selection for initial step');
    end
    if(isnan(t)) keyboard; t=ones(size(t)); end
    if(min(t)==0) keyboard;  t=eps; end;
    if(nargout==0)  keyboard; t=min(t); end
end

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

