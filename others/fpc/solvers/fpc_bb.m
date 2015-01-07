function Out = fpc_bb(n,A,b,mu,M,opts)

% Fixed Point Continuation (FPC) for l1-Regularized Least Squares
%   with BB Steps and Non-monotone Line Search (on full f)
%
% Solves
%
%   min ||x||_1 + (mu/2)*||Ax - b||_M^2.
%
% Inputs
%
% A - an explicit m x n matrix, or a function handle that implements 
%     A*x and A'*x.  If the latter, this function must have the form
%
%       function y = name(x,mode)
%
%     where if mode == 1 then x is n x 1 and y = A*x, and if mode == 2 then 
%     x is m x 1 and y = A'*x.
%
% b - an m x 1 vector.
%
% M - any m x m positive definite matrix, or the empty matrix.  If M
%     is empty, fpc assumes M = I, which reduces the second term of the 
%     objective to (mu/2)*||Ax - b||_2^2 (standard least squares). 
%
% fpc_bb_opts.m describes the available options.  If opts is passed as 
% empty, fpc_bb_opts will be called to obtain the default values.
%
% Outputs
%
%   Out.x    - x at last iteration
%   Out.f    - vector of function values
%   Out.C    - vector of reference function values
%   Out.lam  - vector of ||x||_1
%   Out.step - vector of norm(x - xp)
%   Out.mus  - vector of mu values for each outer iteration
%   Out.itr  - number of iterations to convergence (or Inf if reach mxitr)
%   Out.itrs - vector of number of inner iterations completed during each
%              outer iteration
%   Out.tau   - vector of tau values for each iteration
%   Out.alpha - vector of alpha values for each iteration
%   Out.n2re - if opts.xs exists, is vector of norm(x - xs)/norm(xs).
%              starts with 0th iteration.
%
% Copyright (C) 2007-2008 Elaine Hale, Wotao Yin and Yin Zhang
% FPC is distributed under the GNU GPL, see README.txt and COPYING.txt.
%
% Last modified 2 February 2008

% problem dimension
m = length(b);

% unify implementation of A
if ~isa(A,'function_handle'), A = @(x,mode) explicitMatrix(A,x,mode); end

% same for M
if ~isempty(M) && ~isa(M,'function_handle')
	M = @(x,mode) explicitMatrix(M,x,mode);
end

% get or check opts
opts = fpc_bb_opts(opts);

% check scaling
if opts.scale
    disp('FPC-BB is checking the scaling of your problem because opts.scale = true.');
    [mu,A,b] = scaleProblem(n,mu,A,b,M);
end

% calculate A'*M*b
if isempty(M), AtMb = A(b,2); else AtMb = A(M(b,1),2); end

% check for 0 solution
if mu <= 1/norm(AtMb,'inf'); 
    Out.x = zeros(n,1); Out.itr = 0; Out.itrs = 0; 
    Out.tau = 0; Out.alpha = 0;  Out.mus = mu; Out.lam = 0; Out.step = []; 
    if isempty(M), Out.f = (mu/2)*(b'*b); 
    else Out.f = (mu/2)*(b'*M(b,1)); end
    if isfield(opts,'xs'), Out.n2re = 1; end
    return
end

% initialize x, nu, tau, mu
muf = mu;               % final value of mu
[x,tauD,mu] = fpc_bb_init(n,m,b,AtMb,M,opts);
if mu > muf, mu = muf; end
Out.mus = mu;

% initialize Out.n2re
if isfield(opts,'xs'), xs = opts.xs; else xs = []; end
if ~isempty(xs), Out.n2re = norm(x - xs)/norm(xs); end

gam = opts.gamma;
xtol = opts.xtol;
gtol = opts.gtol;

% prepare for iterations
Out.step = []; Out.itrs = []; Out.f = []; Out.lam = [];
Out.itr = Inf; oitr = 0; 
Out.tau = []; Out.alpha = []; Out.C = []; gp = [];

[g,f,lam,Ax] = get_g(x,mu,A,b,M,AtMb);
Q = 1; C = f;
Out.f = [Out.f; f]; Out.C = [Out.C; C]; Out.lam = [Out.lam; lam];
warned = false;

% main loop
for i = 1:opts.mxitr
    
    % get bb step
    if ~isempty(gp)
        dg = g - gp;
        tau = max((xxp'*xxp)/max(real(xxp'*dg),eps),1);
    else
        tau = tauD;
    end
    nu = tau/mu;    
    
    % store old point
    xp = x; gp = g; Axp = Ax;
    
    % take fixed-point step
    y = x - tau*g; 
    x = sign(y).*max(0,abs(y)-nu);
    
    % calculate g, f, Ax
    [g,f,lam,Ax] = get_g(x,mu,A,b,M,AtMb);
    
    % line search
    alpha = 1; 
    dx = x - xp;
    dg = g - gp;
    dAx = Ax - Axp;
    const = opts.c*((sign(xp) + mu*gp)'*dx);
    if ~isreal(const) && ~warned
        warning('FPC-BB is designed for problems with real data.'); 
        warned = true;
    end
    % xnew = xp + alpha*(x - xp)
    cnt = 0;
    while f > C + alpha*const
        if cnt == 100
            % give up and take safe step
            disp('cnt at 100');
            tau = tauD; alpha = 1; 
            y = xp - tau*gp; 
            x = sign(y).*max(0,abs(y)-nu);
            [g,f,lam,Ax] = get_g(x,mu,A,b,M,AtMb);
            break;
        end
        alpha = alpha*opts.beta;
        [x,g,f,lam,Ax] = update_g(alpha,mu,xp,dx,gp,dg,Axp,dAx,M,b);
        cnt = cnt + 1;
    end
    Qp = Q; Q = gam*Qp + 1; C = (gam*Qp*C + f)/Q;
    Out.f = [Out.f; f]; Out.C = [Out.C; C]; Out.lam = [Out.lam; lam];
    Out.tau = [Out.tau; tau]; Out.alpha = [Out.alpha; alpha];
    
    xxp = x - xp;
    nrmxxp = norm(xxp);
    Out.step = [Out.step; nrmxxp]; 
    
    if ~isempty(xs), Out.n2re = [Out.n2re; norm(x - xs)/norm(xs)]; end
    
    crit1 = nrmxxp/max(norm(xp),1);
    crit2 = mu*norm(g,'inf') - 1;
    
    if (crit1 < xtol*sqrt(muf/mu)) && (crit2 < gtol)
        
        oitr = oitr + 1;
        
        if isempty(Out.itrs), Out.itrs = i;
        else Out.itrs = [Out.itrs; i - sum(Out.itrs)]; end
        
        % stop if reached muf
        if mu == muf
            Out.x = x; Out.itr = i;
            return 
        end
        
        % update mu
        mup = mu;
        if opts.fullMu
            phi = sqrt((2/mu)*(f - lam));
            mu = getNextMu(n,phi,opts.eta,opts.kappa);
        else
            mu = opts.eta*mu;
        end
        mu = min(mu,muf); Out.mus = [Out.mus; mu];
        gp = []; f = (f - lam)*(mu/mup) + lam;
        Q = 1; C = f;
    end
end

% did not converge within opts.mxitr
Out.x = x;
if isempty(Out.itrs), Out.itrs = i;
else Out.itrs = [Out.itrs; i - sum(Out.itrs)]; end

end % fpc

%--------------------------------------------------------------------------
% SUBFUNCTION FOR INITIALIZATION
%--------------------------------------------------------------------------
%
% OUTPUTS -----------------------------------------------------------------
% x    - initialized based on opts.x0 and opts.init.  if opts.x0 exists, 
%       x = opts.x0, otherwise, opts.init determines x:
%           0 - x = zeros(n,1)
%           1 - x = tau*||AtMb||_Inf * ones(n,1)
%           2 - x = tau*AtMb
% tauD - equals opts.tauD if exists, otherwise min(1.999,-1.665*m/n+2.665)
% mu   - set based on x = 0, mu = 1/norm(AtMb,inf), and getNextMu
%--------------------------------------------------------------------------

function [x,tauD,mu] = fpc_bb_init(n,m,b,AtMb,M,opts)

% initialize tau
if isfield(opts,'tauD'), tauD = opts.tau;
else tauD = min(1.999,-1.665*m/n + 2.665); end

% initialize x
if isfield(opts,'x0')
    x = opts.x0;
    if length(x) ~= n, error('User supplied x0 is wrong size.'); end
else
    switch opts.init
        case 0, x = zeros(n,1);
        case 1, x = tauD*norm(AtMb,inf)*ones(n,1);
        case 2, x = tauD*AtMb;
    end
end

% initialize mu
if isfield(opts,'mu0')
    mu = opts.mu0;
elseif opts.fullMu && isfield(opts,'kappa')
    % phi = ||Ax - b||_M
    if isempty(M), phi = norm(b); else phi = sqrt(b'*M(b,1)); end
    mu = getNextMu(n,phi,opts.eta,opts.kappa);
else
    if opts.fullMu
        warning('opts.fullMu = true, but opts.kappa is not supplied.  Switching to opts.fullMu = false.');
        opts.fullMu = false;
    end
    mu = opts.eta/norm(AtMb,inf);
end

end % fpc_bb_init

%--------------------------------------------------------------------------
% SUBFUNCTION FOR CALCULATING NEXT mu
%--------------------------------------------------------------------------
%
% Calculates the next value of mu based on taking a predictor step along
% the pareto curve phi(lam).  The derivative of this curve is derived in 
% 
% van den Berg, E. and M. Friedlander.  In pursuit of a root.  Preprint,
%   2007.
%
% The steplength is chosen so that phi(current)/phi(next) \approx eta.
% Mu is chosen so that the true phi(next) is guaranteed to be <= the 
% predicted value.
%
% INPUTS ------------------------------------------------------------------
% n     - length of x
% phi   - ||Ax - b||_M
% g     - A'M(Ax - b)
% eta   - parameter for choosing how much phi should decrease during next
%         outer iteration.  getNextMu is approximately the same as choosing
%         mu = eta*mu.
% kap   - condition number of A'MA.  See getM_mu.
%--------------------------------------------------------------------------

function mu = getNextMu(n,phi,eta,kap)

mu = eta*sqrt(n*kap)/phi;

end % getNextMu

%--------------------------------------------------------------------------
% SUBFUNCTION FOR CALCULATING g
%--------------------------------------------------------------------------

function [g,f,lam,Ax] = get_g(x,mu,A,b,M,AtMb)

% get A*x
Ax = A(x,1);

% calc g
if isempty(M)
    g = A(Ax,2) - AtMb;
else
    g = A(M(Ax,1),2) - AtMb;
end

% calc f
r = Ax - b; lam = norm(x,1);
if isempty(M)
    f = 0.5*mu*norm(r)^2 + lam;
else
    f = 0.5*mu*r'*M(r,1) + lam;
end

end % get_g

%--------------------------------------------------------------------------
% SUBFUNCTION FOR UPDATING g
%--------------------------------------------------------------------------

function [x,g,f,lam,Ax] = update_g(alpha,mu,xp,dx,gp,dg,Axp,dAx,M,b)

x = xp + alpha*dx;
Ax = Axp + alpha*dAx;

r = Ax - b; lam = norm(x,1);
if isempty(M)
    f = 0.5*mu*norm(r)^2 + lam;
else
    f = 0.5*mu*r'*M(r,1) + lam;
end  
g = gp + alpha*dg;

end % update_g