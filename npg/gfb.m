function out = gfb(F,Gi,C,Lip,xInit,opt)
%   gfb solves a sparse regularized problem using GFB algorithm
%
%                          F(x) + \sum_i G_i(x)
%
%   where F(x) is a function with Lipschitz continuous gradient, G_i(x) are
%   a group of functions with simple proximal operators.  For the input, we
%   require F to be a function handle that can be used in the following
%   form:
%                            [f,grad] = F(x);
%
%   where f and grad are the value and gradient of F at x, respectively.
%
%   Gi is a cell with the i-th element Gi{i} being a structure: Gi{i}.val
%   is a function handle, such that Gi{i}.val(x) returns the value of
%   G_i(x); Gi{i}.prox is another function handle, such that Gi{i}.prox(a)
%   returns the proximal of G_i(x) at a, i.e., the solution of
%
%                    minimize 0.5*||x-a||_2^2 + G_i(x)
%
%   C is a function handle that makes sure that after each iteration
%   C.prox(x) is within the domain of F(x).  C is optional and can be set
%   to C.prox=@(x)x to take no effect.
%
%   xInit       Initial value for estimation of x
%   Lip         Required, is the Lipschitz constant of gradient of F(x)
%
%   opt         The optional structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Example:
%
%     C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
%     opt.debugLevel=2; opt.outLevel=1; opt.m=350;
%     opt.snr=inf; opt.maxItr=1e5;
%     [y,Phi,Phit,Psi,Psit,opt]=loadLinear(opt);
%     opt.u = 10^(-5)*pNorm(Psit(Phit(y)),inf);
%     NLL=@(x) Utils.linearModel(x,Phi,Phit,y);
%     % opt.gamma in [0, 2];
%     opt.gamma=1.9; opt.lambda=1/min(1.5,0.5+1/(opt.gamma/opt.L));
%     G.exact=true;
%     G.val=@(x) opt.u*norm(Psit(x),1);
%     G.prox=@(a,v) Psi(Utils.softThresh(Psit(a),v*opt.u));
%     out=gfb(NLL,{G,C},C,opt.L,zeros(size(opt.trueX)),opt);
%
%   See also: pnpg pds sparseProximal slGaussEx
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
% opt.gamma in [0, 2]
if(~isfield(opt,'gamma')) opt.gamma=1.8; end
% opt.lambda in [0, 1]
if(~isfield(opt,'lambda')) opt.lambda=1; end
if(~isfield(opt,'minItr')) opt.minItr=10; end
if(~isfield(opt,'errorType')) opt.errorType=1; end

if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
% print iterations every opt.verbose lines.
if(~isfield(opt,'verbose')) opt.verbose=1000; end

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
% opt.gamma in [0,2]
% opt.lambda in [0,1]
gamma=opt.gamma/Lip; lambda=min(1.5,0.5+1/gamma)*opt.lambda;

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
    x=C.prox(x);
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
        debug.print(2,sprintf(' %14.12g',cost));
        if(isfield(opt,'trueX'))
            debug.print(2,sprintf(' %12g',RMSE));
        end
        debug.print(2,sprintf(' %12g',difX));
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
out.lambda=lambda; out.w=w; out.gamma=gamma;
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
    fprintf('\t lambda=%g, gamma=%g, w=%g\n',lambda,gamma,w(1));
end

end

