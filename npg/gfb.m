function out = gfb(F,Gi,Lip,xInit,opt)
%   gfb solves a sparse regularized problem using GFB algorithm
%
%                        F(x) + u*r(x) + I_C(x)                     (1)
%
%   where F(x) is the likelihood function, r(x) is the regularization
%   term with u as the regularization parameter, and I_C(x) is an indicator
%   function.   For the input, we require F to be a function handle that
%   can be used in the following form:
%                       [f,grad,hessian] = F(x);
%
%   where f and grad are the value and gradient of F at x, respectively.
%   "hessian" is a function handle that can be use to calculate H*x and
%   x'*H*x with hessian(x,1) and hessian(x,2), respectively.  If Hessian
%   matrix is not available, simply return hessian as empty: [];
%   (See npg/sparseProximal.m and utils/Utils.m for examples)
%
%   The "proximal" parameter is served as a structure with "iterative",
%   "op" and "val" to solve the following subproblem:
%
%                         0.5*||x-a||_2^2+u*r(x)                      (2)
%
%   where proximal.iterative must be false, proximal.prox is exact with no
%   iterations, i.e., the proximal operator has analytical solution:
%                           x=proximal.prox(a,u);
%   where "u" is optional in case r(x) is an indicator function.
%   One example for the case when r(x)=||Psit(x)||_1 with orthogonal Psi,
%   the proximal operator is
%
%           proximal.prox=@(a,t) Psi(softThresh(Psit(a),t));
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
%                       provide Lipschitz constant of F by opt.Lip;
%       debugLevel      An integer value to tune how much debug information
%                       to show during the iterations;
%       outLevel        An integer value to control how many quantities to
%                       put in "out".
%
%   See also: pnpg pds sparseProximal
%
%   Reference:
%   	H. Raguet, J. Fadili, and G. Peyré, “A generalized forward-backward
%       splitting,” SIAM J. Imag. Sci., vol. 6, no. 3, pp. 1199-1226, 2013.
%
%   Author: Renliang Gu (gurenliang@gmail.com)

if(~exist('opt','var') || ~isfield(opt,'u') || isempty(opt.u))
    opt.u=1e-4;
end
if(~exist('Lip','var') || isnan(Lip))
    error('Require Lip as Lipschitz constant of F(x)');
end

% Threshold for relative difference between two consecutive x
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'errorType')) opt.errorType=1; end

if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
% print iterations every opt.verbose lines.
if(~isfield(opt,'verbose')) opt.verbose=100; end

% Output options and debug information
% >=0: minimum output with only results,
% >=1: some cheap output,
% >=2: more detail output and expansive (impairs CPU time, only for debug)
if(~isfield(opt,'outLevel')) opt.outLevel=0; end
if(~isfield(opt,'saveXtrace') || opt.outLevel<2) opt.saveXtrace=false; end

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

% print start information
if(debug.level(2))
    fprintf('\n%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Generalized Forward-Backward (GFB) Method');
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','Itr');
    str=sprintf([str ' %14s'],'Objective');
    if(isfield(opt,'trueX'))
        str=sprintf([str ' %12s'], 'Error');
    end
    str=sprintf([str ' %12s'], '|dx|/|x|');
    str=sprintf([str ' %12s'], '|d Obj/Obj|');
    str=sprintf([str '\t u=%g'],opt.u);
    fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
end

tStart=tic;

% initialization for GFB
itr=0; convThresh=0; x=xInit;
w=zeros(length(Gi),1);

% default to be 1.8/Lipschitz
gamma=1.8/Lip; lambda=1;

[cost,grad]=F(x);
for i=1:length(Gi)
    cost=cost+Gi{i}.val(x);
    z{i}=x;
    w(i)=1/length(Gi);
end

if((opt.outLevel>=1 || debug.level(2)) && isfield(opt,'trueX'))
    RMSE=computError(x);
end

if(opt.outLevel>=1) out.debug={}; end

while(true)

    if(itr >= opt.maxItr || (convThresh>0 && itr>=opt.minItr))
        break;
    end

    itr=itr+1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start of one GFB step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:length(Gi)
        z{i}=z{i}+lambda*(Gi{i}.prox(2*x-z{i}-gamma*grad, gamma/w(i))-x);
    end
    preX=x;
    x=0;
    for i=1:length(Gi)
        x=x+w(i)*z{i};
    end
    preCost=cost;
    [cost, grad]=F(x);
    for i=1:length(Gi)
        cost=cost+Gi{i}.val(x);
    end

    %z1=z1+lambda*(proximal.prox(2*x-z1-gamma*grad,gamma*u/w)-x);
    %z2=z2+lambda*(    obj.prj_C(2*x-z2-gamma*grad      )-x);
    %newX  =w*z1+(1-w)*z2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  end of one GFB step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    difX = relativeDif(x,preX);

    if(difX<=opt.thresh )
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
    end

    if(debug.level(2))
        debug.print(2,sprintf(' %5d',itr));
        debug.print(2,sprintf(' %14.8g',cost));
        if(isfield(opt,'trueX'))
            debug.print(2,sprintf(' %12g',RMSE));
        end
        debug.print(2,sprintf(' %12g',difX));
        debug.print(2,sprintf(' %12g', difCost));
        debug.clear_print(2);
        if(mod(itr,opt.verbose)==0) debug.println(2); end
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
out.lambda=lambda; out.w=w; out.gamma=gamma;
if(opt.outLevel>=2)
    out.grad=grad;
end
if(debug.level(1))
    fprintf('\nCPU Time: %g, objective=%g',toc(tStart),cost);
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

end

