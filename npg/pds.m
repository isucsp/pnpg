function out = pds(F,H,Psi,Psit,G,Lip,xInit,opt)
%   pds solves a sparse regularized problem using PDS algorithm
%
%                      F(x) + r(Psit(x)) + I_C(x)                     (1)
%
%   where F(x) is the likelihood function, r(s) is the regularization
%   term including u as the regularization parameter, and I_C(x) is an
%   indicator function.   For the input, we require F to be a function
%   handle that can be used in the following form:
%
%                       [f,grad,hessian] = F(x);
%
%   where f and grad are the value and gradient of F at x, respectively.
%   "hessian" is a function handle that can be use to calculate H*x and
%   x'*H*x with hessian(x,1) and hessian(x,2), respectively.  If Hessian
%   matrix is not available, simply return hessian as empty: [];
%   (See npg/sparseProximal.m and utils/Utils.m for examples)
%
%   The "H" parameter is served as a structure with "exact",
%   "op" and "val" to solve the following subproblem:
%
%                         0.5*||x-a||_2^2+u*r(x)                      (2)
%
%   where H.exact must be false ??? not true , H.prox is exact with no
%   iterations, i.e., the H operator has analytical solution:
%                           x=H.prox(a,u);
%   where "u" is optional in case r(x) is an indicator function.
%   One example for the case when r(x)=||Psit(x)||_1 with orthogonal Psi,
%   the H operator is
%
%           H.prox=@(a,t) Psi(softThresh(Psit(a),t));
%
%   H.val(x) returns the value of r(x).
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
%    1. L. Condat, “A primal-dual splitting method for convex optimization
%       involving Lipschitzian, proximable and linear composite terms,” J.
%       Optim. Theory Appl. , vol. 158, no. 2, pp. 460-479, 2013.
%    2. B. C. Vũ, “A splitting algorithm for dual monotone inclusions
%       involving cocoercive operators,” Adv. Comput. Math., vol. 38, no.
%       3, pp. 667-681, 2013.
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

% In case of projection as H
if(nargin(H.proxConj)==1)
    proximalOp=H.proxConj;
    H.proxConj=@(a,u) proximalOp(a);
end

% print start information
if(debug.level(2))
    fprintf('\n%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Primal-Dual Splitting (PDS) Method');
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

% initialization for PDS
itr=0; convThresh=0; x=xInit; y=zeros(size(x));

sigma=(sqrt(Lip^2/16+1)-Lip/4)*1;
tau=1/(sigma+Lip/2) * 1;
rho=(2-Lip/2/(1/tau-sigma)) * 1;

[f,grad]=F(x);
Psit_x=Psit(x);
cost=f+G.val(x)+H.val(Psit_x);

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
    %  start of one PDS step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ybar=H.proxConj(y+sigma*Psit_x,sigma);
    xbar=G.prox(x-tau*(grad+Psi(2*ybar-y)),tau);
    preX=x;
    x=x+rho*(xbar-x);
    y=y+rho*(ybar-y);
    Psit_x=Psit(x);

    [f,grad]=F(x);
    preCost=cost;
    cost = f+G.val(x)+H.val(Psit_x);

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  end of one PDS step  %
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
out.sigma=sigma; out.tau=tau; out.rho=rho;
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

