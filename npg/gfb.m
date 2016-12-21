function out = gfb(NLL,proximal,xInit,opt)
%   gfb solves a sparse regularized problem using GFB algorithm
%
%                        NLL(x) + u*r(x) + I_C(x)                     (1)
%
%   where NLL(x) is the likelihood function, r(x) is the regularization
%   term with u as the regularization parameter, and I_C(x) is an indicator
%   function.   For the input, we require NLL to be a function handle that
%   can be used in the following form:
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
%   where proximal.iterative must be false, proximal.op is exact with no
%   iterations, i.e., the proximal operator has analytical solution:
%                           x=proximal.op(a,u);
%   where "u" is optional in case r(x) is an indicator function.
%   One example for the case when r(x)=||Psit(x)||_1 with orthogonal Psi,
%   the proximal operator is
%
%           proximal.op=@(a,t) Psi(softThresh(Psit(a),t));
%
%   proximal.penalty(x) returns the value of r(x).
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
%   See also: pnpg pds sparseProximal
%
%   Reference:
%   	H. Raguet, J. Fadili, and G. Peyré, “A generalized forward-backward
%       splitting,” SIAM J. Imag. Sci., vol. 6, no. 3, pp. 1199-1226, 2013.
%
%   Author: Renliang Gu (gurenliang@gmail.com)

% default to not use any constraints.
if(~exist('opt','var') || ~isfield(opt,'prj_C') || isempty(opt.prj_C))
    opt.prj_C=@(x) x;
end

if(~isfield(opt,'u')) opt.u=1e-4; end
if(~isfield(opt,'Lip')) opt.Lip=nan; end

if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'saveXtrace')) opt.saveXtrace=false; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive x
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'errorType')) opt.errorType=1; end

% Output options and debug information
% 0: minimum output, 1: some output, 2: more detail output
if(~isfield(opt,'outLevel')) opt.outLevel=0; end
if(~isfield(opt,'NLL_Pen') || opt.outLevel<=1) opt.NLL_Pen=false; end

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
if(opt.outLevel>=4)
    figCost=1000; figure(figCost);
end
if(opt.outLevel>=5)
    figRes=1001; figure(figRes);
end
if(opt.outLevel>=6)
    figX=1002; figure(figX);
end

if(proximal.iterative)
    error('proximal.op need to be an analytical solution!');
end

% In case of projection as proximal
if(nargin(proximal.op)==1)
    proximalOp=proximal.op;
    proximal.op=@(a,u) proximalOp(a);
end

% print start information
if(debug.level(1))
    fprintf('%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Generalized Forward-Backward (GFB) Method');
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

tStart=tic;

itr=0; convThresh=0; x=xInit; difX=1;
cost=NLL(x)+opt.u*proximal.penalty(x);

z1=x; z2=x;

% default to be 1.8/Lipschitz
r=1.8/opt.Lip;
lambda=1; w=0.5;

if(opt.outLevel>=1) out.debug={}; end

while(true)
    itr=itr+1;
    %if(mod(itr,100)==1 && itr>100) save('snapshotFST.mat'); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start of one GFB step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [oldCost,grad] = NLL(x);
    z1=z1+lambda*(proximal.op(2*x-z1-r*grad,r*u/w)-x);
    z2=z2+lambda*(obj.prj_C(2*x-z2-r*grad)-x);
    newX  =w*z1+(1-w)*z2;

    NLLVal=NLL(newX);
    penVal = proximal.penalty(newX);
    cost = NLLVal+opt.u*penVal;

    difX = relativeDif(x,newX);
    x = newX;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  end of one GFB step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    out.cost(itr) = cost;
    if(itr>1)
        difCost=abs(cost-out.cost(itr-1))/cost;
    end
    if(opt.outLevel>=1)
        out.time(itr)=toc(tStart);
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
        out.stepSize(itr) = lambda;
        if(opt.NLL_Pen)
            out.NLLVal(itr)=NLLVal;
            out.penVal(itr)=penVal;
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

end

