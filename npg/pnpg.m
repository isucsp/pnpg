function out = pnpg(NLL,proximal,xInit,opt)
%   pnpg solves a sparse regularized problem
%
%                           NLL(x) + u*r(x)                           (1)
%
%   where NLL(x) is the likelihood function and r(x) is the regularization
%   term with u as the regularization parameter.   For the input, we
%   require NLL to be a function handle that can be used in the following
%   form:
%                       [f,grad,hessian] = NLL(x);
%
%   where f and grad are the value and gradient of NLL at x, respectively.
%   "hessian" is a function handle that can be use to calculate H*x and
%   x'*H*x with hessian(x,1) and hessian(x,2), respectively.  If Hessian
%   matrix is not available, simply return hessian as empty: [];
%   (See npg/sparseProximal.m and utils/Utils.m for examples)
%
%   The "proximal" parameter is served as a structure with "iterative",
%   "op" and "penalty" to solve the following subproblem:
%
%                         0.5*||x-a||_2^2+u*r(x)                      (2)
%
%   1. When proximal.iterative=true, proximal.op is iterative and should be
%   called in the form of
%             [x,itr,p]=proximal.op(a,u,thresh,maxItr,pInit);
%   where
%       pInit           initial value of internal variable (e.g., dual
%                       variable), can be [] if not sure what to give;
%       maxItr          maximum number of iteration to run;
%       x               solution of eq(2);
%       itr             actual number of iteration run when x is returned;
%       p               value of internal variable when terminates, can be
%                       used as the initial value (pInit) for the next run.
%
%   2. When proximal.iterative=false, proximal.op is exact with no
%   iterations, i.e., the proximal operator has analytical solution:
%                           x=proximal.op(a,u);
%   where "u" is optional in case r(x) is an indicator function.
%
%   proximal.penalty(x) returns the value of r(x).
%   (See npg/sparseProximal.m for an example)
%
%   xInit       Initial value for estimation of x
%
%   opt         The optional structure for the configuration of this algorithm (refer to
%               the code for detail)
%       u               Regularization parameter;
%       initStep        Method for the initial step size, can be one of
%                       "hessian", "bb", and "fixed".  When "fixed" is set,
%                       provide Lipschitz constant of NLL by opt.Lip;
%       debugLevel      An integer value to tune how much debug information
%                       to show during the iterations;
%       outLevel        An integer value to control how many quantities to
%                       put in "out".
%
%   Reference:
%   Author: Renliang Gu (gurenliang@gmail.com)

% default to not use any constraints.
if(~exist('opt','var') || ~isfield(opt,'prj_C') || isempty(opt.prj_C))
    opt.prj_C=@(x) x;
end

if(~isfield(opt,'u')) opt.u=1e-4; end

if(~isfield(opt,'stepIncre')) opt.stepIncre=0.9; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'initStep')) opt.initStep='hessian'; end
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'saveXtrace')) opt.saveXtrace=false; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive x
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'Lip')) opt.Lip=nan; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'errorType')) opt.errorType=1; end
if(~isfield(opt,'gamma')) gamma=2; else gamma=opt.gamma; end
if(~isfield(opt,'b')) b=0.25; else b=opt.b; end

if(~isfield(opt,'relInnerThresh')) opt.relInnerThresh=1e-2; end
if(~isfield(opt,'cumuTol')) opt.cumuTol=4; end
if(~isfield(opt,'incCumuTol')) opt.incCumuTol=true; end
if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=true; end
if(~isfield(opt,'maxInnerItr')) opt.maxInnerItr=100; end
if(~isfield(opt,'maxPossibleInnerItr')) opt.maxPossibleInnerItr=1e3; end

% Output options and debug information
% 0: minimum output, 1: some output, 2: more detail output
if(~isfield(opt,'outLevel')) opt.outLevel=0; end
if(~isfield(opt,'NLL_Pen') || opt.outLevel<=1) opt.NLL_Pen=false; end
if(~isfield(opt,'collectTheta') || opt.outLevel<=1) opt.collectTheta=false; end

if(isfield(opt,'trueX'))
    switch opt.errorType
        case 0
            trueX = opt.trueX/pNorm(opt.trueX);
            computError= @(xxx) 1-(realInnerProd(xxx,trueX)^2)/sqrNorm(xxx);
        case 1
            trueXNorm=sqrNorm(opt.trueX);
            if(trueXNorm==0) trueXNorm=eps; end
            computError = @(xxx) sqrNorm(xxx-opt.trueX)/trueXNorm;
        case 2
            trueXNorm=pNorm(opt.trueX);
            if(trueXNorm==0) trueXNorm=eps; end
            computError = @(xxx) pNorm(xxx-opt.trueX)/trueXNorm;
    end
end

debug=Debug(opt.debugLevel);
if(opt.outLevel>=4)
    figCost=1000; figure(figCost);
end
if(opt.outLevel>=5)
    figRes=1001; figure(figRes);
end
if(opt.outLevel>=6)
    figX=1002; figure(figX);
end

if(proximal.iterative && abs(nargout(proximal.op))>1 && opt.outLevel>0)
    collectInnerItr=true;
else
    collectInnerItr=false;
end

% In case of projection as proximal
if(nargin(proximal.op)==1)
    proximalOp=proximal.op;
    proximal.op=@(a,u) proximalOp(a);
end

% print start information
if(debug.level(1))
    fprintf('%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Projected Nestrov''s Proximal-Gradient (PNPG) Method');
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','Itr');
    str=sprintf([str ' %14s'],'Objective');
    if(isfield(opt,'trueX'))
        str=sprintf([str ' %12s'], 'Error');
    end
    str=sprintf([str ' %12s %4s'], '|dx|/|x|', 'αSrh');
    str=sprintf([str ' %12s'], '|d Obj/Obj|');
    str=sprintf([str '\t u=%g'],opt.u);
    fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
end

t=stepSizeInit(opt.initStep,opt.Lip);

tic;

itr=0; convThresh=0; theta=1;
preT=t; preX=0; x=xInit; difX=1;
cost=NLL(x)+opt.u*proximal.penalty(x);

if(opt.outLevel>=1) out.debug={}; end
if(proximal.iterative) pInit=[]; end
if(opt.adaptiveStep) cumu=0; end

while(true)
    itr=itr+1;
    %if(mod(itr,100)==1 && itr>100) save('snapshotFST.mat'); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start of one PNPG step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    numLineSearch=0; goodMM=true;
    if(opt.adaptiveStep)
        incStep=false;
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
        [oldCost,grad] = NLL(xbar);
    end

    % start of line Search
    while(true)
        numLineSearch = numLineSearch+1;

        if(opt.adaptiveStep)
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
        end

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
            if(numLineSearch<=20 && t>0)
                t=t/opt.stepShrnk;
                % Penalize if there is a step size increase just now
                if(opt.adaptiveStep && incStep)
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
    newPen = proximal.penalty(newX);
    newObj = newCost+opt.u*newPen;

    if(newObj>cost)
        if(goodMM)
            if(restart()) continue; end
            if(runMore()) continue; end
        else
            reset(); % both theta and step size;
            if(runMore()) continue; end
        end
        % give up and force it to converge
        debug.appendLog('_ForceConverge');
        innerItr=0;
        preX=x; difX=0;
    else
        if(proximal.iterative)
            pInit=pInit_;
            innerItr=innerItr_;
        else
            innerItr=0;
        end
        difX = relativeDif(x,newX);
        preX = x;
        x = newX;
        theta = newTheta;
        cost = newObj;
        NLLVal=newCost;
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

    out.cost(itr) = cost;
    if(itr>1)
        difCost=abs(cost-out.cost(itr-1))/cost;
    end
    if(opt.outLevel>=1)
        out.time(itr)=toc;
        out.difX(itr)=difX;
        if(itr>1) out.difCost(itr)=difCost; end
        if(~isempty(debug.log()))
            out.debug{size(out.debug,1)+1,1}=itr;
            out.debug{size(out.debug,1),2}=debug.log();
            debug.clearLog();
        end;
    end

    if(opt.outLevel>=2)
        out.numLineSearch(itr) = numLineSearch;
        out.stepSize(itr) = 1/t;
        if(opt.NLL_Pen)
            out.NLLVal(itr)=NLLVal;
            out.penVal(itr)=penVal;
        end
        if(opt.collectTheta) out.theta(itr)=theta; end
        if(collectInnerItr) out.innerItr(itr)=innerItr; end;
        if(debug.level(2))
            out.BB(itr,1)=stepSizeInit('BB');
            out.BB(itr,2)=stepSizeInit('hessian');
        end
        if(opt.saveXtrace) out.xTrace(:,itr)=x; end
    end

    debug.print(1,sprintf(' %5d',itr));
    debug.print(1,sprintf(' %14.12g',cost));
    if(isfield(opt,'trueX'))
        out.RMSE(itr)=computError(x);
        debug.print(1,sprintf(' %12g',out.RMSE(itr)));
    end
    debug.print(1,sprintf(' %12g %4d',difX,numLineSearch));
    if(itr>1)
        debug.print(1,sprintf(' %12g', difCost));
    else
        debug.print(1,sprintf(' %12s', ' '));
    end
    debug.clear_print(1);
    if(mod(itr,opt.verbose)==0) debug.println(1); end

    if(itr>1 && opt.outLevel>=4)
        set(0,'CurrentFigure',figCost);
        if(isfield(opt,'trueX')) subplot(2,1,1); end
        if(out.cost(itr)>0)
            semilogy(itr-1:itr,out.cost(itr-1:itr),'k'); hold on;
            title(sprintf('cost(%d)=%g',itr,out.cost(itr)));
        end

        if(isfield(opt,'trueX'))
            subplot(2,1,2);
            semilogy(itr-1:itr,out.RMSE(itr-1:itr)); hold on;
            title(sprintf('RMSE(%d)=%g',itr,out.RMSE(itr)));
        end
        drawnow;
    end

    if(opt.NLL_Pen && itr>1 && opt.outLevel>=5)
        set(0,'CurrentFigure',figRes);
        subplot(2,1,1);
        semilogy(itr-1:itr,out.NLLVal(itr-1:itr),'r'); hold on;
        subplot(2,1,2);
        semilogy(itr-1:itr,out.penVal(itr-1:itr),'b'); hold on;
        drawnow;
    end

    if(opt.outLevel>=6)
        set(0,'CurrentFigure',figX); showImgMask(x,opt.mask);
        drawnow;
    end

    if(itr>1 && difX<=opt.thresh )
        convThresh=convThresh+1;
    end

    if(itr >= opt.maxItr || (convThresh>2 && itr>opt.minItr))
        break;
    end
end
out.x=x; out.itr=itr; out.opt = opt; out.date=datestr(now);
if(opt.outLevel>=2)
    out.grad=grad;
end
if(debug.level(0))
    fprintf('\nCPU Time: %g, objective=%g',out.time(end),out.cost(end));
    if(isfield(opt,'trueX'))
        fprintf(', RMSE=%g\n',out.RMSE(end));
    else
        fprintf('\n');
    end
end

function reset()
    theta=1; t=min([t;max(stepSizeInit('hessian'))]);
    debug.appendLog('_Reset');
    debug.printWithoutDel(1,'\t reset');
end
function res=restart()
    % if has monmentum term, restart
    res=pNorm(xbar-x,1)~=0;
    if(res)
        theta=1;
        debug.appendLog('_Restart');
        debug.printWithoutDel(1,'\t restart');
        itr=itr-1;
    end
end
function res=runMore()
    res=false;
    if(~proximal.iterative) return; end
    if(innerItr_<opt.maxInnerItr && opt.relInnerThresh>1e-6)
        opt.relInnerThresh=opt.relInnerThresh/10;
        debug.printWithoutDel(1,...
            sprintf('\n decrease relInnerThresh to %g',...
            opt.relInnerThresh));
        itr=itr-1;
        res=true;
        %% IMPORTANT! if not requir too high accuracy
        %% use 1e3 for maxInnerItr
    elseif(innerItr_>=opt.maxInnerItr &&...
            opt.maxInnerItr<opt.maxPossibleInnerItr)
        opt.maxInnerItr=opt.maxInnerItr*10;
        debug.printWithoutDel(1,...
            sprintf('\n increase maxInnerItr to %g',opt.maxInnerItr));
        itr=itr-1;
        res=true;
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

    % if(~isempty(dfx) && abs(fx-fy)/max(max(fx,fy),1) < 1e-10)
    %     % In this case, use stronger condition to avoid numerical issue
    %     test=(realInnerProd(x_minus_y,dfx-dfy) <= L*sqrNorm(x_minus_y)/2);
    % else
        test=((fx-fy)<=realInnerProd(x_minus_y,dfy)+L*sqrNorm(x_minus_y)/2);
    % end
end

