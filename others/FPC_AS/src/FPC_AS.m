function [x, Out] = FPC_AS(n,A,b,mu,M,opts)

%--------------------------------------------------------------------------
% GENERAL DESCRIPTION & INPUTS
%-------------------------------------------------------------------------- 
%
% [x,Out] = fpc_AS(n,A,b,mu,M,opts)
%
% Solves
%
%  (P) min f(x) = mu*||x||_1 + 0.5*||Ax - b||_M^2.
%
% n the dimension of x, the column rank of A
%
% A either an explicit m x n matrix or an A_operator object representing
%   a matrix implicitly. When A is relatively small and/or operations A*x 
%   and A.'*x cannot be computed much faster through certain means, it is 
%   recommend that A is provided explicitly. Otherwise, A should be created
%   as an A_operator object, the implementation of which is provided with 
%   the solver. 
%
%   To create an A_operator object, two functions or function handles for
%   computing A*x and A.'*x, respectively, must be given. Suppose they
%   are AA and AT,
%       AA: a function handle such that AA(x) = A*x,
%       AT: a function handle such that AT(x) = A.*x.
%   Then A can be created as an A_operator by
%       A = A_operator(@(x) AA(x), @(x) AT(x));
%
%   An example for A being an implicit matrix that performs a discrete 
%   cosine transform (DCT) and returns the subset of the DCT result 
%   corresponding to omega is
%
%     function y=pdct(x,picks); y=dct(x); y=y(picks); end
%     function y=pidct(x,n,picks); y=zeros(n,1); y(picks)=x; y=idct(y); end
%     A = A_operator(@(x) pdct(x,omega), @(x) pidct(x,n,omega));  
%
% b must be an m x 1 vector, the row rank of A must equal to m
%
% M either an m x m positive definite matrix or the empty matrix [].  
%   If M=[], FPCAS treats M = I, which reduces the last term in (P) 
%   to 0.5*||Ax - b||_2^2. 
%   
%   FPCAS usually works better if the maximum eigenvalue of A' M A is 
%   close to, but no larger than, 1. If it is larger than 1, the fact that 
%   0.5*||Ax - b||_M^2  will dominate ||x||_1  will affect the performance
%   of the algorithm. To remove this effect, either scale the input argument
%   A by theta*A for some appropriate theta<1 or let FPCAS do so by 
%   setting opts.scale_A = 1 below.
%
% opts a structure of options. It is an optional argument, so it can be 
%   ignored or set empty.  Some of the frequently used fields  include:
%   
%   'mxitr': max number of iterations
%           default: 1000, valid range: [1, 100000]
%   'gtol': termination criterion on ``crit2'', the maximum norm of sub-gradient
%           default: 1e-08, valid range: [0, 1]
%   'gtol_scale_x': termination criterion on ``crit2'' scaled by max(norm(x), 1). 
%           default: 1e-12, vaild range: [0, 1]
%   'f_value_tol': Tolerance on the optimal objective value, stop if f(x) <= f_value_tol
%           default: 0, vaild range: [0, inf]
%   'sub_mxitr': max number of iterations in each sub-optimization.
%       	default: 80, valid range: [1, 100000]
%   'sub_opt_meth': choice of sub-optimization methods.
%           default: 'lbfgs', valid values: {'lbfgs','lbfgsb','pcg'};
%   'lbfgs_m': storage number of L-BFGS. 
%           default: 5, valid range: [1, 100]
%   'scale_A': on/off switch for scaling the input matrix A so that the max of eigs(A*A.') equals 1. 
%           default: 0, valid range: {0, 1}
%   'minK': an estimate of the number of nonzero components in optimal solution. 
%           default: m/2, valid range: [1, n]
%   'zero': a lower bound of minimal magnitude of nonzero components of optimal solution
%           default: 1e-08, valid range: [0, 1e+10]
%   'dynamic_zero': on/off switch for setting `zero' dynamically
%           default: 0, valid range: {0, 1}
%   'xs': optimal solution for non-algorithmic purposes such as progress display 
%           default: empty vector, valid range: [-Inf, Inf]
%   'record': print level, -1=quiet, 0=some output, 1=more output. 
%           default: 0 valid range: {-1, 0, 1}
%   'PrintOptions': print options, 0=quiet, 1=output 
%           default: 0, valid range: {0, 1}            
%
%   To set a field of opts, do  opts.[fieldname] = value. 
%   A complete list of options are given in the manual.
%
%   *****************************IMPORTANT*********************************
%   The default values of these options are configured to solve (P) to a
%   high accuracy. If 'b' is contaminated by noise or if it is NOT
%   necessary to solve (P) very accurately, please try, for example, 
%        opts.sub_mxitr = 10; opts.gtol = 1e-3; opts.gtol_scale_x = 1e-6;
%   ***********************************************************************
%
%--------------------------------------------------------------------------
% OUTPUTS
%--------------------------------------------------------------------------
% x             -exit solution
% Out           -a structure having the following fields
%    .cpu       -the total amount of CPU time 
%    .exit      -exit status, 1='normal', 10='exceed max iterations'
%    .mesg      -a string describe the detailed exit status
%    .itr       -the number of iterations
%    .f         -exit function value, mu*||x||_1 + 0.5*||Ax - b||_M^2
%    .nrm1x     -exit l1 norm ||x||_1
%    .rNorm     -exit l2 discrepancy term ||Ax - b||_M
%    .g         -exit gradient of 0.5*||Ax - b||_M^2
%
%    .zero      -final thresholding value for computing termination criteria 
%    .crit2     -violation of optimality, norm(g(T)-mu,'inf'), where T = union(nz_x, z_xa),
%                nz_x = x>Out.zero, the set of nonzero components whose magnitude of x 
%                       are larger than Out.zero
%                z_xa = ~nz_x & (||g|-mu| > Out.zero); the set of nonzero components 
%                       whose magnitude of |g|-mu are larger than Out.zero
%
%    .nnz_x     -number of elements such that x_i > Out.zero
%    .nCont     -number of continuation steps taken
%    .nSubOpt   -number of sub-optimization problems solved
%    .nProdA    -total number of A*x    
%    .nProdAt   -total number of A'*x
%    .nfe       -number of A*x called by shrinkage
%    .nge       -number of A'*x called by shrinkage
%    .nfe_sub   -number of A*x called by sub-optimization
%    .nge_sub   -number of A'*x called by sub-optimization
%
%    .opts      -options used
%
%    --> The following fields are available if an optimal solution opts.xs 
%    is provided in the input. opts.xs is compared to x after x is truncated
%    in the way that x_i=0 if x_i <Out.zero. The comparison results are given in
%
%    .sgn       -number of nonzero components of x that have different sign compared to those of opts.xs,
%    .miss      -number of (missed) components that are zero in x but nonzero in opts.xs
%    .over      -number of (overshot) components that are nonzero in x but zero in opts.xs
%
%    They are computed as 
%         jnt  = union(nz_xs, nz_x); nz_xs = abs(xs) > opts.eps;
%         sgn  = nnz(sign(x(jnt))~=sign(xs(jnt))); 
%         miss = nnz(nz_xs&(~nz_x))
%         over = nnz(nz_x&(~nz_xs))
% 
%--------------------------------------------------------------------------
% Copyright (c) 2008.  Zaiwen Wen and Wotao Yin
%
% version 1.1, 10/2008, Zaiwen Wen
%--------------------------------------------------------------------------

    function y = AAt(x)
        if isempty(M); y = A*((A')*x); else y = A*(M*((A')*x)); end; nAATEigs = nAATEigs + 1;
    end
    function scale_mat_A()
        nAATEigs = 0; smax = 1;  opts_eig.tol = opts.tol_eig; opts_eig.disp = 0; opts_eig.issym = true; opts_eig.maxit = opts.eigs_mxitr;

        err_eig.identifier = []; try; smax = sqrt(eigs(@AAt,m,1,'lm',opts_eig)); catch err_eig = lasterror;  end;
        if strcmp(err_eig.identifier, 'MATLAB:eigs:complexFunction') == 1
            opts_eig.isreal = 0; try; smax = sqrt(eigs(@AAt,m,1,'lm',opts_eig)); catch err_eig = lasterror;  end;
        end
        if smax ~= 1

            A = A*(1/smax); b = (1/smax)*b; opts.tauD = max(opts.tauD, 2/smax);
            if opts.record == 1; fprintf('The max eigen of A*A^T is %3.2e, tauD: %3.2e\n', smax, opts.tauD); end
        end
    end

    function val = calc_f()
        if isempty(M); qx = 0.5*norm(res)^2; else qx = 0.5*res'*(M*res); end;  val = mu*nrm1x + qx;
    end

    function grad = calc_g()
        if isempty(M); grad = A'*res;  else  grad = (A')*(M*res); end; Out.nge = Out.nge+1;
    end

    function Ax_x_update()
        Ax = A*x; res=Ax-b; norm_x = max(norm(x),1); nrm1x = norm(x,1); Out.nfe = Out.nfe+1;
    end

    function Ax_x_update_ls()
        Ax = Axp + alpha*Ad; res=Ax-b; if isempty(M); qx = 0.5*norm(res)^2; else qx = 0.5*res'*(M*res); end
    end

    function x_prev_save()
        xp = x; gp = g; Axp = Ax; resp = res; fp = f; qxp = qx; nrm1xp = nrm1x;
    end

    function compute_crits()
        crit1 = nrmxxp/norm_x; Out.rNorm = sqrt(qx*2);



        nz_x = abs(x)> epsIx; nnz_x = nnz(nz_x);   z_x = ~nz_x;


        abs_g_mu = abs(g)/mu-1;  z_x_a = z_x & (abs_g_mu > epsIx);
        nrm_g_mu = norm(abs_g_mu(nz_x), inf); crit2 = mu*max([ nrm_g_mu, norm(abs_g_mu(z_x_a), inf)]);
        f_rel = abs( (fp-f)/max(abs([f,fp,1])));

        y = x - g; nrm_d1x = norm(sign(y).*max(0,abs(y)-mu)-x);  crittol = 1/kappa_g_d * nrm_d1x;


        if ~opts.dynamic_zero;  epsIx = opts.zero; else

            rIx = mu*norm(x.*abs_g_mu) + min(crit2, crit2f);

            epsIx = min([max([0.001*sqrt(rIx)/norm_x, opts.zero]), nrm1x/n]);


        end


        if muf == mu; crit2f = crit2; else
            abs_g_muf = abs(g)/muf-1;  z_x_a = z_x & (abs_g_muf > epsIx);
            nrm_g_muf = norm(abs_g_muf(nz_x), inf);     crit2f = muf*max([ nrm_g_muf, norm(abs_g_muf(z_x_a), inf)]);
        end

        if (opts.record==1);
            fprintf('itr=%d mu=%3.2e crit(2,2f)=(%3.2e,%3.2e), grad/nrmxxp=%3.2e, res=%3.2e,',itr,mu,crit2/norm_x,crit2f/norm_x,norm(alpha*tau*g(nz_x))/(nrmxxp), Out.rNorm/norm_x);
            if ~isempty(xs),
                jnt=nz_xs&nz_x;
                fprintf(' err=%4.2e nnz=%d sgn=%d miss=%d over=%d\n',norm(x-xs)/norm(xs),nnz_x,nnz(sign(x(jnt))~=sign(xs(jnt))),nnz(nz_xs&(~nz_x)),nnz(nz_x&(~nz_xs)));
            else fprintf(' nnz=%d\n',nnz_x);  end
            Ix_Diff_g = (sign(x) == sign(g)); nIx_Diff_g = nnz(Ix_Diff_g(nz_x));
            crit_compl = norm(x.*abs_g_mu);
            fprintf('sgn_x_g: %d, rIx: %3.2e, epsIx: %4.3e, comp: %3.2e, nrmx: %3.2e, crittol: %3.2e, kappa: %3.2e\n\n', nIx_Diff_g, rIx, epsIx, crit_compl/norm_x, norm_x, crittol, kappa_g_d);
        end
    end

    function final_solu_info()
        Out.cpu = cputime - Out.cpu; Out.nfe_sub = nfe_sub; Out.nge_sub = nge_sub; Out.nProdA = Out.nProdA + nfe_sub + Out.nfe; Out.nProdAt = Out.nProdAt + nge_sub + Out.nge;
        Out.nrm1x = nrm1x; Out.f = muf*nrm1x + qx; Out.g = g; Out.itr = itr; Out.opts = opts; Out.crit2 = crit2f; Out.nnz_x = nnz_x;
        if opts.record >= 0
            fprintf(1,'\n-----------------------------------------------------\n');
            fprintf(1,'solver: %s\n', mfilename('fullpath'));  if opts.PrintOptions; PrintOption(opts); end
            if ~isempty(xs),
                jnt=nz_xs&nz_x;

                Out.ab_err =  norm(x-xs,inf); Out.rel_err = norm(x-xs)/norm(xs); Out.sgn = nnz(sign(x(jnt))~=sign(xs(jnt))); Out.miss =  nnz(nz_xs&(~nz_x)); Out.over = nnz(nz_x&(~nz_xs));
                fprintf('Problem information: n=%d, m=%d, K=%d\n',n,m,nnz(xs));
                fprintf('abs.err=%4.2e, rel.err=%4.2e, nnz(x)=%d, sgn=%d, miss=%d, over=%d\n', Out.ab_err, Out.rel_err, Out.nnz_x, Out.sgn, Out.miss, Out.over);
            else fprintf('Problem information: n=%d, m=%d, nnz(x)=%d\n',n,m,nnz(x)); end
            fprintf('time: %6.3fs, crit2: %3.2e, nrm1x: %3.2e, ||Ax-b||: %3.2e, f: %3.2e\n', Out.cpu, Out.crit2, nrm1x, Out.rNorm, Out.f);
            fprintf('cost: num. of shrinkage: %d, num. of sub-opt: %d, num. of continuation: %d\n', itr, Out.nSubOpt, Out.nCont);
            fprintf('num. of A*x from (total, shrink, sub-opt): (%d, %d, %d), num. of A''*x: (%d, %d, %d)\n',  Out.nfe + Out.nfe_sub, Out.nfe, Out.nfe_sub,  Out.nge + Out.nge_sub, Out.nge, Out.nge_sub);
            fprintf('Message: %s\n',Out.mesg);
        end
        clear global MT AT cT bb nfe_sub nge_sub;
    end

    function result = fpc_stagnate()
        result = 0;
        if nnz_x < 1 ; result = 0; return; end
        if itr-itr_prev_subopt<opts.min_itr_shrink; result=0; return; end

        if isequal(nz_x,nz_x_prev_subopt) && itr-itr_prev_subopt<opts.max_itr_shrink;  result=0; return; end

        if (alpha*tau*norm(g(nz_x))/nrmxxp>kappa_g_d) && crit2/norm_x < min(opts.tol_start_sub,crittol); result = 1; kappa_g_d = kappa_g_d*kappa_rho; end


        if f_rel < opts.f_rel_tol; result = 1; end;
        if  itr-itr_prev_subopt>opts.max_itr_shrink; result = 1; end;

        if result > 0 && opts.record==1; fprintf('stagnation %d detected.\n', result); end;
    end

    function subopt()
        cT = muf*sign(x(nz_x))'; nfe_sub_p = nfe_sub;
        if ~isempty(M); MT = M(nz_x,nz_x); end
        if isa(A,'A_operator'); AT.A = A; AT.col = nz_x; AT.n = n; else  AT = A(:,nz_x); end

        x(z_x)=0; gp_good=false; nz_x_prev_subopt = nz_x; itr_prev_subopt = itr;
        if strcmp(opts.sub_opt_meth, 'lbfgsb') == 1
            lb = repmat(-inf,nnz_x,1); lb(cT>0) = 0;
            ub = repmat(inf,nnz_x,1);  ub(cT<0) = 0;

            lbfgsb_callback = []; factrtol = min(1e-6,0.1*opts.gtol); pgtol = min(1e-6,0.1*opts.gtol);





            nb = repmat(2,nnz_x,1);  nb(cT>0) = 1;  nb(cT<0) = 3;
            optsl = lbfgs_options('iprint', -1, 'maxits', opts.sub_mxitr ,'m',opts.lbfgs_m, 'factr', factrtol, 'pgtol',pgtol);
            [x(nz_x), f, gsub, Outsub]= flbfgsb(x(nz_x), cT, AT, b, MT, lb, ub, nb, optsl); nfe_sub = nfe_sub+ Outsub.nfe; nge_sub = nge_sub+ Outsub.nge;
            res = gsub.res; Ax = gsub.Ax; qx = gsub.qx; norm_x = max(norm(x),1);  nrm1x = norm(x,1); f = mu*nrm1x + qx; if isobject(A); g = gsub.g; else g = calc_g(); end
        elseif strcmp(opts.sub_opt_meth, 'lbfgs') == 1
            optsl.debug = false; optsl.m = opts.lbfgs_m; optsl.maxitr = opts.sub_mxitr; optsl.gtol = min(1e-6,0.1*opts.gtol); optsl.xtol = min(1e-6,0.1*opts.gtol);
            [x(nz_x), f, gsub, Outsub] = fminLBFGS_Q(x(nz_x), cT, AT, bb, MT, optsl); nfe_sub = nfe_sub+ Outsub.nfe; nge_sub = nge_sub+ Outsub.nge;

            res = gsub.res; Ax = gsub.Ax; qx = gsub.qx; norm_x = max(norm(x),1);  nrm1x = norm(x,1); f = mu*nrm1x + qx; if isobject(A); g = gsub.g; else g = calc_g(); end
        elseif strcmp(opts.sub_opt_meth, 'pcg') == 1
            if isstruct(AT);
                if isempty(MT); bT = (AT.A')*bb;  else bT = (AT.A')*(MT*bb);  end; cT = bT(AT.col)- cT';
            else
                if isempty(MT);  cT = AT'*bb - cT';   else cT = AT'*(MT*bb) - cT'; end
            end
            nge_sub = nge_sub+1;
            optsl.maxitr = opts.sub_mxitr; optsl.gtol = min(1e-6,0.1*opts.gtol);
            [x(nz_x), Outsub.msg,Outsub.relres,Outsub.nfe] = wpcg(@AsubCG, cT,optsl.gtol,optsl.maxitr,[],[],x(nz_x));

            Ax_x_update(); f = calc_f();  g = calc_g();
        end
        if (opts.record==1); fprintf('sub-opt: nfe: %3d, fp: %3.2e, f: %3.2e, f-rel: %3.2e, fp-f: %3.2e, g*d: %3.2e\n', nfe_sub - nfe_sub_p, fp, f, abs( (fp-f)/max(abs([f,fp,1]))), fp - f, gp'* (x-xp)); end
    end


Out.cpu = cputime;


if (nargin < 4), error('At least 4 parameters are required! SYNOPSIS: [x, Out] = FPC_AS(n,A,b,mu,M,opts)');  end
if (nargin < 5), M = []; opts=[]; end
if (nargin < 6), opts = [];  end
if mu < 0; error('mu = %e < 0. Please choose mu > 0', mu); end
if ~isa(A,'numeric') && ~isa(A,'A_operator'); error('A must be a numerical matrix or A_operator\n'); end

muf = mu;  m = length(b); opts = FPC_AS_opts(opts,n,m);
nAATEigs = 0; if opts.scale_A;  scale_mat_A(); end; Out.nProdA = nAATEigs; Out.nProdAt = nAATEigs;


if isempty(M), AtMb = (A')*b; else  AtMb = (A')*(M*b); end; normAtMb = norm(AtMb,'inf'); Out.nProdAt = Out.nProdAt + 1;

xs = []; if isfield(opts,'xs'), if ~isempty(opts.xs); xs = opts.xs; nz_xs = abs(xs)>opts.eps; end; end
gamma = opts.gamma;


global MT AT cT bb nfe_sub nge_sub;
bb = b; if isempty(M); MT=[]; end

x = zeros(n,1);  xp = zeros(n,1); g = zeros(n,1);  gp = zeros(n,1); res = zeros(m,1);  qx = 0; Ax = []; Axp = [];
nz_x = false(n,1); z_x = false(n,1); nnz_x = 0;
crit1 = inf; crit2 = inf; crit2f = inf; norm_x = inf; nrm1x = inf; nrmxxp = inf; nrm_d1x = inf; crittol = inf;
Out.nfe = 0; Out.nge = 0; itr = 0; alpha = 0; f = inf; fp = inf; f_rel = inf;
rIx = inf; epsIx = opts.zero; Out.nCont = 0;  Out.nSubOpt = 0; nfe_sub = 0; nge_sub = 0;

itr_prev_subopt = 0; itr_prev_mu_dec = 0; nz_x_prev_subopt = false(n,1);
kappa_g_d = opts.kappa_g_d; kappa_rho = opts.kappa_rho;

if mu >= normAtMb; 
    x = zeros(n,1);     Out.crit2 = 0; Out.exit = 1; Out.mesg = 'zero is optimal'; 
    compute_crits();    final_solu_info();      return; 
end

tauD = opts.tauD; tau = tauD; 
if ~isempty(opts.x0);
    x = opts.x0;
else
    switch opts.init;  case 0, x = zeros(n,1); case 1, x = normAtMb*ones(n,1); case 2, x = AtMb; end
end

Ax_x_update(); g=calc_g(); compute_crits();

if nnz_x <= m && crit2f < opts.gtol ;
    Out.exit = 1; Out.mesg = 'initial x optimal'; final_solu_info();  return;
else
    mu_prev = normAtMb;    mu = max(mu_prev*opts.eta,10*muf);
end
f=calc_f(); fp = f;


gp_good = false; Q = 1; C = f;
for itr = 1:opts.mxitr

    if (gp_good); tau = (xxp'*xxp)/(xxp'*(g - gp)); tau = max(min(tau, opts.tau_max), opts.tau_min); else tau = tauD; end


    x_prev_save(); gp_good=true;


    y = x - tau*g;  x = sign(y).*max(0,abs(y)-tau*mu); Ax_x_update();   f = calc_f();


    alpha = 1; dx = x - xp;


    const = opts.c*(gp'*dx + mu.*(nrm1x - nrm1xp));

    if strcmp(opts.ls_meth, 'nmls') == 1
        cnt = 0;
        dols = false; if f > C + alpha*const;  dols = true; Ad = Ax - Axp; c1 = Ad'*Ad; c2 = Ad'*resp; end
        while dols
            if cnt == opts.ls_mxitr

                tau = tauD;
                y = xp - tau*gp; x = sign(y).*max(0,abs(y)-tau*mu); Ax_x_update(); f = calc_f();
                break;
            end
            alpha = alpha*opts.beta; tau = tau*opts.beta;

            x = xp + alpha*dx; nrm1x = norm(x,1); f = 0.5*alpha^2*c1 + alpha*c2 + qxp + mu*nrm1x;
            if f <= C + alpha*const; dols = false; Ax_x_update_ls(); break; end;
            cnt = cnt + 1;
        end
    elseif strcmp(opts.ls_meth, 'Hybridls') == 1

        if f > C + alpha*const
            Ad = Ax - Axp; c1 = Ad'*Ad; c2 = Ad'*resp;  parsl1.debug = false;

            [alpha, f1, x, Outl1] = SolveOneDimL1(n, 0, 1, xp, dx, mu, c1, c2, parsl1);

            if abs(alpha-1) > 1e-16
                f = f1 + qxp; nrm1x = norm(x,1); Ax_x_update_ls();
            end
        end
    end
    g = calc_g();

    Qp = Q; Q = gamma*Qp + 1; C = (gamma*Qp*C + f)/Q;

    if (opts.record==1); fprintf('tau: %3.2e, alpha: %3.2e, fp: %3.2e, f: %3.2e, f-rel: %3.2e, fp-f: %3.2e, qx/nrm1x: %3.2e\n', tau, alpha, fp, f, abs( (fp-f)/max(abs([f,fp,1]))), fp - f, qx/nrm1x); end


    xxp = x - xp; nrmxxp = norm(xxp);  compute_crits();

    if (crit2f<opts.gtol) || (crit2f/norm_x<opts.gtol_scale_x) || (mu == muf && f_rel < opts.f_rel_tol);
        Out.exit = 1; Out.mesg = 'shrinkage optimal'; final_solu_info();  return;
    elseif f < opts.f_value_tol
        Out.exit = 5; Out.mesg = 'the current function value is less than opts.f_value'; final_solu_info();  return;
    end


    fpc_stall = fpc_stagnate();
    if fpc_stall == 1
        if nnz_x > m;
            xnzx = sort(abs(x),'descend'); epsIx = xnzx(opts.minK);
            if opts.record==1; fprintf('check sparsity: nnz_x: %6d, zero: %4.3e, |x|(m/2): %4.3e\n', nnz_x, epsIx, xnzx(opts.minK)); end
            nz_x = abs(x)>epsIx; nnz_x = nnz(nz_x); z_x = ~nz_x; x(z_x) = 0; Ax_x_update(); f = calc_f(); g = calc_g();
        end;
        x_prev_save();


        subopt(); Out.nSubOpt = Out.nSubOpt + 1;  compute_crits();
        if (crit2f<opts.gtol)  || (crit2f/norm_x<opts.gtol_scale_x)  || (mu == muf && f_rel < opts.f_rel_tol)
            Out.exit = 1; Out.mesg = 'sub-opt optimal'; final_solu_info();  return;
        elseif f < opts.f_value_tol
            Out.exit = 5; Out.mesg = 'the current function value is less than opts.f_value'; final_solu_info();  return;
        end
    end


    if ( (crit2/norm_x<opts.gtol_scale_x) || (crit2f < crit2 ) || fpc_stall > 0 ) && (mu > muf)

        mu_prev = mu; kappa_g_d = opts.kappa_g_d; max_g_zx = max(abs(g(z_x)));
        if (itr - itr_prev_mu_dec <= opts.min_itr_shrink)
            mu = mu_prev *opts.eta;
        end
        mu = max([min([max_g_zx*opts.eta,mu*opts.eta]), muf]);
        itr_prev_mu_dec = itr;
        Out.nCont = Out.nCont + 1; f = mu*nrm1x + qx; gp_good=false; Q=1; C=f; fp = f;
        if opts.record==1; fprintf('Cont %d-th: muf: %4.3e, mu_prev=%4.3e, mu=%4.3e, eta: %4.3e\n', Out.nCont, muf, mu_prev, mu, opts.eta); end
    end

end

Out.exit = 10; Out.mesg = 'exceed max iterations'; final_solu_info();
end
