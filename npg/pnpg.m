function out = pnpg(NLL,proximal,xInit,opt)
%   pnpg solves a sparse regularized problem
%
%                           NLL(x) + u*r(x)                           (1)
%
%   where x is a real vector value to optimize over, NLL(x) is the
%   likelihood function, and r(x) is the regularization term with u as the
%   regularization parameter.   For the input, we require NLL to be a
%   function handle that can be used in the following form:
%
%                       [f,grad,hessian] = NLL(x);
%
%   where f and grad are the value and gradient of NLL at x, respectively.
%   "hessian" is a function handle that can be use to calculate H*x and
%   x'*H*x with hessian(x,1) and hessian(x,2), respectively.  If Hessian
%   matrix is not available, simply return hessian as empty: []; (See
%   npg/sparseProximal.m and utils/Utils.m for examples)
%
%   The "proximal" parameter is served as a structure with "exact",
%   "op" and "val" to solve the following subproblem:
%
%                         0.5*||x-a||_2^2+u*r(x)                      (2)
%
%   1. When proximal.exact=false, proximal.prox is inexact and should be
%   called in the form of
%             [x,itr,p]=proximal.prox(a,u,thresh,maxItr,pInit);
%   where
%       pInit           initial value of internal variable (e.g., dual
%                       variable), can be [] if not sure what to give;
%       maxItr          maximum number of iteration to run;
%       x               solution of eq(2);
%       itr             actual number of iteration run when x is returned;
%       p               value of internal variable when terminates, can be
%                       used as the initial value (pInit) for the next run.
%
%   2. When proximal.exact=true, proximal.prox is exact with no
%   iterations, i.e., the proximal operator has analytical ??? solution:
%                           x=proximal.prox(a,u);
%   where "u" is optional in case r(x) is an indicator function.
%
%   proximal.val(x) returns the value of r(x).
%   (See npg/sparseProximal.m for an example)
%
%   xInit       Initial value for estimation of x
%
%   opt         The optional structure for the configuration of this algorithm (refer to
%               the code for detail)
%       prj_C           A function handle to project signal x to a convex set
%                       C;
%       u               Regularization parameter;
%       initStep        Method for the initial step size, can be one of
%                       "hessian", "bb", and "fixed".  When "fixed" is set,
%                       provide Lipschitz constant of NLL by opt.Lip;
%       debugLevel      An integer value to tune how much debug information
%                       to show during the iterations;
%       outLevel        An integer value to control how many quantities to
%                       put in "out".
%
%   Usage:
%
%
%   Reference:
%   Author: Renliang Gu (gurenliang@gmail.com)

% default to not use any constraints.
if(~exist('opt','var') || ~isfield(opt,'prj_C') || isempty(opt.prj_C))
    opt.prj_C=@(x) x;
end

if(~isfield(opt,'u')) opt.u=1e-4; end

if(~isfield(opt,'stepIncre')) opt.stepIncre=0.5^0.2; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
% By default disabled.  Remember to use a value around 5 for the Poisson model
% with poor initialization.
if(~isfield(opt,'preSteps')) opt.preSteps=0; end
if(~isfield(opt,'initStep')) opt.initStep='hessian'; end
% Threshold for relative difference between two consecutive x
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'Lip')) opt.Lip=inf; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'maxLineSearch')) opt.maxLineSearch=20; end
if(~isfield(opt,'errorType')) opt.errorType=1; end
if(~isfield(opt,'gamma')) gamma=2; else gamma=opt.gamma; end
if(~isfield(opt,'b')) b=0.25; else b=opt.b; end

if(~isfield(opt,'dualGap')) opt.dualGap=false; end
if(~isfield(opt,'relInnerThresh')) opt.relInnerThresh=1e-2; end
% decend rate of (epsion*theta)^2 is O(1/i^epsilonDecRate)
if(~isfield(opt,'epsilonDecRate')) opt.epsilonDecRate=1; end
if(~isfield(opt,'cumuTol')) opt.cumuTol=4; end
if(~isfield(opt,'incCumuTol')) opt.incCumuTol=true; end
if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=true; end
if(~isfield(opt,'maxInnerItr')) opt.maxInnerItr=100; end
if(~isfield(opt,'maxPossibleInnerItr')) opt.maxPossibleInnerItr=1e3; end
if(strcmpi(opt.initStep,'fixed') && isinf(opt.Lip))
    error('The "fixed" for opt.initStep can not be used together with infinite opt.Lip');
end

% Debug output information
% >=0: no print,
% >=1: only report results,
% >=2: detail output report, 
% >=4: plot real time cost and RMSE,
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
% print iterations every opt.verbose lines.
if(~isfield(opt,'verbose')) opt.verbose=100; end

% Output options and debug information
% >=0: minimum output with only results,
% >=1: some cheap output,
% >=2: more detail output and expansive (impairs CPU time, only for debug)
if(~isfield(opt,'outLevel')) opt.outLevel=0; end
if(~isfield(opt,'saveXtrace') || opt.outLevel<2) opt.saveXtrace=false; end
if(~isfield(opt,'collectOtherStepSize') || opt.outLevel<2)
    opt.collectOtherStepSize=false;
end

if(isfield(opt,'trueX'))
    switch opt.errorType
        case 0
            trueX = opt.trueX/pNorm(opt.trueX,2);
            computError= @(xxx) 1-(innerProd(xxx,trueX)^2)/sqrNorm(xxx);
        case 1
            trueXNorm=sqrNorm(opt.trueX);
            if(trueXNorm==0) trueXNorm=eps; end
            computError = @(xxx) sqrNorm(xxx-opt.trueX)/trueXNorm;
        case 2
            trueXNorm=pNorm(opt.trueX,2);
            if(trueXNorm==0) trueXNorm=eps; end
            computError = @(xxx) pNorm(xxx-opt.trueX,2)/trueXNorm;
    end
end

debug=Debug(opt.debugLevel);
if(debug.level(4))
    figCostRMSE=1000; figure(figCostRMSE);
end

% In case of projection as proximal
if(nargin(proximal.prox)==1)
    proximalOp=proximal.prox;
    proximal.prox=@(a,u) proximalOp(a);
end

% print start information
if(debug.level(2))
    fprintf('\n%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Projected Nestrov''s Proximal-Gradient (PNPG) Method');
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','Itr');
    str=sprintf([str ' %14s'],'Objective');
    if(isfield(opt,'trueX'))
        str=sprintf([str ' %12s'], 'Error');
    end
    str=sprintf([str ' %12s %4s'], '|dx|/|x|', 'lSrh');
    if(~proximal.exact)
        str=sprintf([str ' %4s'], 'iItr');
    end
    str=sprintf([str ' %12s'], '|d Obj/Obj|');
    str=sprintf([str '\t u=%g'],opt.u);
    fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
end

tStart=tic;

itr=0; convThresh=0; x=xInit; theta=1; preX=x; itrRes=0;
NLLVal=NLL(x);
penVal=proximal.val(x);
cost=NLLVal+opt.u*penVal;
goodStep=true;
t=stepSizeInit(opt.initStep,opt.Lip);

if((opt.outLevel>=1 || debug.level(2)) && isfield(opt,'trueX'))
    RMSE=computError(x);
end

if(opt.outLevel>=1) out.debug={}; end
if(~proximal.exact)
    pInit=[];
    difX=1;
end
if(opt.adaptiveStep) cumu=0; end

while(true)
    if(itr>=opt.maxItr || convThresh>=3) break; end
    itr=itr+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start of one PNPG step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numLineSearch=0; goodMM=true; incStep=false;
    if(opt.adaptiveStep)
        if(cumu>=opt.cumuTol)
            % adaptively increase the step size
            t=t*opt.stepIncre;
            cumu=0;
            incStep=true;
        end
    else
        %newTheta=(1+sqrt(1+4*B*theta^2))/2;
        if(itr==1)
            newTheta=1;
        else
            newTheta=1/gamma+sqrt(b+theta^2);
        end
        xbar=x+(theta-1)/newTheta*(x-preX);
        xbar=opt.prj_C(xbar);
        [oldNLL,grad] = NLL(xbar);
    end

    % start of line Search
    while(true)
        numLineSearch = numLineSearch+1;

        if(opt.adaptiveStep)
            if(itr==1)
                newTheta=1;
            else
                B=t/preT;
                newTheta=1/gamma+sqrt(b+B*theta^2);
                %newTheta=(1+sqrt(1+4*B*theta^2))/2;
            end
            xbar=x+(theta-1)/newTheta*(x-preX);
            xbar=opt.prj_C(xbar);
            [oldNLL,grad] = NLL(xbar);
        end

        if(proximal.exact)
            newX=proximal.prox(xbar-grad/t,opt.u/t);
        else
            if(opt.dualGap)
                [newX,innerItr_,pInit_]=proximal.prox(xbar-grad/t,opt.u/t,...
                    opt.relInnerThresh/((itr-itrRes)^opt.epsilonDecRate)/(newTheta^2),...
                    opt.maxInnerItr,pInit);
            else
                [newX,innerItr_,pInit_]=proximal.prox(xbar-grad/t,opt.u/t,...
                    opt.relInnerThresh*difX,...
                    opt.maxInnerItr,pInit);
            end
        end

        newNLL=NLL(newX);
        if((newNLL-oldNLL)<=innerProd(newX-xbar,grad)+t*sqrNorm(newX-xbar)/2)
            if(itr<=opt.preSteps && opt.adaptiveStep && goodStep)
                cumu=opt.cumuTol;
            end
            break;
        else
            if(numLineSearch<=opt.maxLineSearch && t<opt.Lip)
                t=t/opt.stepShrnk; goodStep=false;
                % Penalize if there is a step size increase just now
                if(incStep)
                    incStep=false;
                    if(opt.incCumuTol)
                        opt.cumuTol=opt.cumuTol+4;
                    end
                end
            else  % don't know what to do, mark on debug and break
                goodMM=false;
                debug.appendLog('_FalseMM');
                %debug.appendLog(sprintf('\nnewCost-oldNLL%10.10g',newNLL-oldNLL));
                %debug.appendLog(sprintf('\ndif=%10.10g',...
                %   newNLL-oldNLL-innerProd(newX-xbar,grad)+t*sqrNorm(newX-xbar)/2));
                break;
            end
        end
    end
    newPen = proximal.val(newX);
    newCost = newNLL+opt.u*newPen;

    % using eps reduces numerical issue around the point of convergence
    if((newCost-cost)>1e-14*norm([newCost,cost],inf))
        if(goodMM)
            if(restart())
                if(opt.dualGap)
                    opt.relInnerThresh=...
                        opt.relInnerThresh/((itr-itrRes)^opt.epsilonDecRate)/(newTheta^2);
                end
                itr=itr-1;
                itrRes=itr;
                continue;
            end
            if(runMore()) itr=itr-1; continue; end
        end
        %reset(); % both theta and step size;

        % give up and force it to converge
        debug.appendLog('_ForceConverge');
        %debug.appendLog(sprintf('\nNewObj-cost=%10.10g',newCost-cost));
        innerItr=0;
        preX=x; difX=0;
        preCost=cost;
    else
        if(~proximal.exact)
            pInit=pInit_;
            innerItr=innerItr_;
        end
        difX = relativeDif(x,newX);
        preX = x;
        x = newX;
        theta = newTheta;
        preCost=cost;
        cost = newCost;
        NLLVal=newNLL;
        penVal=newPen;
    end

    if(opt.adaptiveStep)
        preT=t;
        if(numLineSearch==1)
            cumu=cumu+1;
        else
            cumu=0;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  end of one PNPG step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    if(difX<=opt.thresh && itr>=opt.minItr)
        convThresh=convThresh+1;
    end

    % skip the rest if not needed
    if(opt.outLevel<1 && opt.debugLevel<2)
        continue;
    end

    difCost=abs(cost-preCost)/max(1,abs(cost));
    if(isfield(opt,'trueX'))
        preRMSE=RMSE; RMSE=computError(x);
    end

    if(opt.outLevel>=1)
        out.time(itr)=toc(tStart);
        out.cost(itr)=cost;
        out.difX(itr)=difX;
        out.difCost(itr)=difCost;
        out.theta(itr)=theta;
        out.numLineSearch(itr) = numLineSearch;
        out.stepSize(itr) = 1/t;
        out.NLLVal(itr)=NLLVal;
        out.penVal(itr)=penVal;
        if(~proximal.exact)
            out.innerItr(itr)=innerItr;
        end;
        if(isfield(opt,'trueX'))
            out.RMSE(itr)=RMSE;
        end
        if(~isempty(debug.log()))
            out.debug{size(out.debug,1)+1,1}=itr;
            out.debug{size(out.debug,1),2}=debug.log();
            debug.clearLog();
        end;
    end

    if(opt.outLevel>=2)
        if(opt.saveXtrace) out.xTrace{itr}=x; end
        if(opt.collectOtherStepSize)
            out.BB(itr,1)=stepSizeInit('BB');
            out.BB(itr,2)=stepSizeInit('hessian');
        end
    end

    if(debug.level(2))
        debug.print(2,sprintf(' %5d',itr));
        debug.print(2,sprintf(' %14.8g',cost));
        if(isfield(opt,'trueX'))
            debug.print(2,sprintf(' %12g',RMSE));
        end
        debug.print(2,sprintf(' %12g %4d',difX,numLineSearch));
        if(~proximal.exact)
            debug.print(2,sprintf(' %4d',innerItr));
        end
        debug.print(2,sprintf(' %12g', difCost));
        if(~debug.clear_print(2))
            if(mod(itr,opt.verbose)==0) debug.println(2); end
        end
    end

    if(debug.level(4))
        set(0,'CurrentFigure',figCostRMSE);
        if(isfield(opt,'trueX')) subplot(2,1,1); end
        if(cost>0)
            semilogy(itr-1:itr,[preCost,cost],'k'); hold on;
            title(sprintf('cost(%d)=%g',itr,cost));
        end

        if(isfield(opt,'trueX'))
            subplot(2,1,2);
            semilogy(itr-1:itr,[preRMSE, RMSE]); hold on;
            title(sprintf('RMSE(%d)=%g',itr,RMSE));
        end
        drawnow;
    end
end
out.x=x; out.itr=itr; out.opt = opt; out.date=datestr(now);
if(opt.outLevel>=2)
    out.grad=grad;
end
if(debug.level(1))
    fprintf('\n%s: CPU Time: %g, objective=%g',mfilename,toc(tStart),cost);
    if(isfield(opt,'trueX'))
        if(debug.level(2))
            fprintf(', RMSE=%g\n',RMSE);
        else
            fprintf(', RMSE=%g\n',computError(x));
        end
    else
        fprintf('\n');
    end
end

function reset()
    theta=1;
    t=min([t;stepSizeInit(opt.initStep,opt.Lip)]);
    debug.appendLog('_Reset');
    debug.printWithoutDel(2,'\t reset');
end
function res=restart()
    % if has monmentum term, restart
    res=((pNorm(xbar-x,0)~=0) && ...
        (oldNLL+opt.u*proximal.val(xbar)>newCost));
    if(~res) return; end

    theta=1;
    debug.appendLog('_Restart');
    debug.printWithoutDel(2,'\t restart');
end
function res=runMore()
    res=false;
    if(proximal.exact) return; end
    if(innerItr_<opt.maxInnerItr && opt.relInnerThresh>1e-9)
        opt.relInnerThresh=opt.relInnerThresh/10;
        debug.printWithoutDel(2,...
            sprintf('\n decrease relInnerThresh to %g',...
            opt.relInnerThresh));
        debug.appendLog('_DecInnerThresh');
        res=true;
    elseif(innerItr_>=opt.maxInnerItr &&...
            opt.maxInnerItr<opt.maxPossibleInnerItr)
        opt.maxInnerItr=opt.maxInnerItr*10;
        debug.printWithoutDel(2,...
            sprintf('\n increase maxInnerItr to %g',opt.maxInnerItr));
        debug.appendLog('_IncInnerMaxItr');
        res=true;
    end
end

function t_=stepSizeInit(select,Lip,delta)
    switch (lower(select))
        case 'bb'   % use BB method to guess the initial stepSize
            if(~exist('delta','var')) delta=1e-5; end
            [~,grad1] = NLL(x);
            temp = delta*grad1/pNorm(grad1,2);
            temp = x-opt.prj_C(x-temp);
            [~,grad2] = NLL(x-temp);
            t_ = abs(innerProd(grad1-grad2,temp))/sqrNorm(temp);
        case 'hessian'
            [~,grad1,hessian] = NLL(x);
            if(isempty(hessian))
                if(~exist('delta','var')) delta=1e-5; end
                temp = delta*grad1/pNorm(grad1,2);
                temp = x-opt.prj_C(x-temp);
                [~,grad2] = NLL(x-temp);
                t_ = abs(innerProd(grad1-grad2,temp))/sqrNorm(temp);
            else
                t_ = hessian(grad1,2)/sqrNorm(grad1);
            end
        case 'fixed'
            t_ = Lip;
        otherwise
            error('unkown selection for initial step');
    end
    if(isnan(t_) || t_<=0)
        error('PNPG is having a negative or NaN step size, do nothing and return!!');
    end
end

end

