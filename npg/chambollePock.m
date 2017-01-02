function out = chambollePock(F,G,K,xInit,opt)
% Solves the minimization of the following problem:
%         F(K(x))+G(x)
% whose primal-dual is
%    min_x  max_y  y'*K(x)-F^*(y)+G(x)
% where F^*(x) is the conjugate of F(x).
%
% Note that F.opConj(a,u) solves 0.5*||x-a||_2^2+u*F^*(x)
% F.penalty(x) returns F(x);
% Note that G.op(a,u) solves 0.5*||x-a||_2^2+u*G(x)
% G.penalty(x) returns G(x);
%
% Reference:
%  [1] A. Chambolle and T. Pock, "A first-order primal-dual algorithm for
%      convex problems with applications to imaging,” J. Math. Imaging Vis.,
%      vol. 40, no. 1, pp. 120–145, 2011.
%
% Author: Renliang Gu (gurenliang@gmail.com)
%

% default to not use any constraints.
if(~exist('opt','var') || ~isfield(opt,'prj_C') || isempty(opt.prj_C))
    opt.prj_C=@(x) x;
end

if(~isfield(opt,'u')) opt.u=1e-4; end

if(~isfield(opt,'stepIncre')) opt.stepIncre=0.9; end
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

if(~isfield(opt,'relInnerThresh')) opt.relInnerThresh=1e-2; end
if(~isfield(opt,'cumuTol')) opt.cumuTol=4; end
if(~isfield(opt,'incCumuTol')) opt.incCumuTol=true; end
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

% In case of projection as G
if(nargin(G.op)==1)
    proximalOp=G.op;
    G.op=@(a,u) proximalOp(a);
end

% print start information
if(debug.level(2))
    fprintf('\n%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Chambolle-Pock''s primal-dual Method');
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','Itr');
    str=sprintf([str ' %14s'],'Objective');
    if(isfield(opt,'trueX'))
        str=sprintf([str ' %12s'], 'Error');
    end
    str=sprintf([str ' %12s'], '|dx|/|x|');
    if(G.iterative)
        str=sprintf([str ' %4s'], 'iItr');
    end
    str=sprintf([str ' %12s'], '|d Obj/Obj|');
    str=sprintf([str '\t u=%g'],opt.u);
    fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
end

tStart=tic;

itr=0; convThresh=0; x=xInit;
Kx=K.forward(x); preKx=Kx; Kxbar=Kx;

cost=F.penalty(Kx)+G.penalty(x);

if((opt.outLevel>=1 || debug.level(2)) && isfield(opt,'trueX'))
    RMSE=computError(x);
end

if(opt.outLevel>=1) out.debug={}; end
if(G.iterative)
    pInit=[];
    difX=1;
end

if(F.iterative)
    error('Iterative F.op is not supported yet!');
end

tau=opt.tau;
sigma=opt.sigma;

y=Kx;
if(iscell(Kx))
    for i=1:length(Kx)
        y{i}=sigma*Kx{i};
    end
    y=F.opConj(y,sigma);
else
    y=F.opConj(sigma*Kx,sigma);
end

while(true)

    if(itr >= opt.maxItr || (convThresh>2 && itr>=opt.minItr))
        break;
    end

    itr=itr+1;
    %if(mod(itr,100)==1 && itr>100) save('snapshotFST.mat'); end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  start of one Chambolle-Pock step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % xbar=2*x-preX;
    if(iscell(Kx))
        for i=1:length(Kx)
            Kxbar{i}=2*Kx{i}-preKx{i};
            y{i}=y{i}+sigma*Kxbar{i};
        end
        y=F.opConj(y,sigma);
    else
        Kxbar=2*Kx-preKx;
        y=F.opConj(y+sigma*Kxbar,sigma);
    end

    preX = x;
    if(G.iterative)
        [x,innerItr_,pInit_]=G.op(x-tau*K.backward(y),tau,...
            opt.relInnerThresh*difX*0,                    opt.maxInnerItr,pInit);
    else
        x=G.op(x-tau*K.backward(y),tau);
    end
    preKx = Kx;
    Kx=K.forward(x);

    preCost=cost;
    cost=F.penalty(Kx)+G.penalty(x);
    difX = relativeDif(x,preX);

    if(G.iterative)
        pInit=pInit_;
        innerItr=innerItr_;
    else
        innerItr=0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %  end of one PNPG step  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    if(itr>1 && difX<=opt.thresh )
        convThresh=convThresh+1;
    end

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
        if(G.iterative)
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
        if(opt.saveXtrace) out.xTrace(:,itr)=x; end
        if(opt.collectOtherStepSize)
            out.BB(itr,1)=stepSizeInit('BB');
            out.BB(itr,2)=stepSizeInit('hessian');
        end
    end

    if(debug.level(2))
        debug.print(2,sprintf(' %5d',itr));
        debug.print(2,sprintf(' %14g',cost));
        if(isfield(opt,'trueX'))
            debug.print(2,sprintf(' %12g',RMSE));
        end
        debug.print(2,sprintf(' %12g',difX));
        if(G.iterative)
            debug.print(2,sprintf(' %4d',innerItr));
        end
        debug.print(2,sprintf(' %12g', difCost));
        debug.clear_print(2);
        if(mod(itr,opt.verbose)==0) debug.println(2); end
    end

    if(debug.level(4) && itr>1)
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

end % function out = chambollePock(F,G,K,xInit,opt)

