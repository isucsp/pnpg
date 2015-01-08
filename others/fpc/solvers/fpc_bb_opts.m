function opts = fpc_bb_opts(opts)

% Options for Fixed Point Continuation with BB Steps (FPC-BB)
%
% If opts is empty upon input, opts will be returned containing the default
% options for fpc_bb.m.  
%
% Alternatively, if opts is passed with some fields already defined, those
% fields will be checked for errors, and the remaining fields will be added
% and initialized to their default values.
%
% Table of Options.  ** indicates default value.
%
% FIELD   OPTIONAL  DESCRIPTION
% .x0       YES     Initial value of x.  If not defined, x will be
%                   initialized according to .init.
% .mu0      YES     Initial value of mu.  If not defined, mu will be
%                   initialized by applying standard procedure to 
%                   1/norm(AtMb,inf).
% .xs       YES     Original signal xs.  If passed, fpc will calculate and
%                   output ||x - xs||/||xs|| in vector Out.n2re.
% .init     YES     If .x0 is not defined, .init specifies how x is to be
%                   initialized.  
%                       0 -> zeros(n,1)
%                       1 -> x = tau*||AtMb||_Inf * ones(n,1)
%                    ** 2 -> x = tau*AtMb **
% .xtol     NO      Tolerance on norm(x - xp)/norm(xp).
%                   ** 1E-4 **
% .gtol     NO      Tolerance on mu*norm(g,'inf') - 1
%                   ** 0.2 **
% .tauD      YES    0 < .tauD < 2.  If not specified, tauD is initialized 
%                   using a piecewise linear function of delta = m/n.
% .mxitr    NO      Maximum number of inner iterations.
%                   ** 1000 **
% .eta      NO      Ratio of current ||b - Ax|| to approx. optimal 
%                   ||b - Ax|| for next mu value.
%                   ** 4 **
% .fullMu   NO      If true, then mu = eta*sqrt(n*kap)/phi, where phi = 
%                   ||Ax - b||_M, which guarantees that phi(now)/phi(next) 
%                   >= eta.  Otherwise mu = eta*mu.
%                   ** false **
% .kappa    YES     Required if fullMu.  Is ratio of max and min 
%                   eigenvalues of M^{1/2}AA'M^{1/2} (before scaling).
%                   ** not supplied **
% .gamma    NO      Weighting for reference function value, between 0 and
%                   1.  0 is monotone line search, 1 is average function 
%                   value.
%                   ** 0.85 **
% .c        NO      Constant for Armijo condition.
%                   ** 1E-3 **
% .beta     NO      Step-size reduction factor for Armijo condition.
%                   ** 0.5 **
% .scale    NO      If true, fpc_bb will calculate the maximum eigenvalue  
%                   of A'*M*A and scale mu, A and b appropriately if it is 
%                   not equal to 1.
%                   ** true **
%
% Copyright (C) 2007-2008 Elaine Hale, Wotao Yin, and Yin Zhang
% FPC is distributed under the GNU GPL, see README.txt and COPYING.txt.
%
% Last modified 2 February 2008.

if isfield(opts,'x0')
    if ~isvector(opts.x0) || ~min(isfinite(opts.x0))
        error('If used, opts.x0 should be an n x 1 vector of finite numbers.');
    end
elseif isfield(opts,'init')
    if ~isinInterval(opts.init,0,2,true) || opts.init ~= floor(opts.init)
        error('opts.init must be an integer between 0 and 2, inclusive.');
    end
else
    opts.init = 2; 
end

if isfield(opts,'mu0')
    if ~isinInterval(opts.mu0,0,Inf,false)
        error('If used, opts.mu0 must be a positive scalar.');
    end
end

if isfield(opts,'xs')
     if ~isvector(opts.xs) || ~min(isfinite(opts.xs))
        error('If passed, opts.xs should be an n x 1 vector of finite numbers.');
     end
end

if isfield(opts,'xtol')
    if ~isinInterval(opts.xtol,0,1,false)
        error('opts.xtol is tolerance on norm(x - xp)/norm(xp).  Should be in (0,1).');
    end
else
    opts.xtol = 1E-4;
end

if isfield(opts,'gtol')
    if ~isinInterval(opts.gtol,0,1,false)
        error('opts.gtol is tolerance on mu*norm(g,''inf'') - 1.  Should be in (0,1).');
    end
else
    opts.gtol = 0.2;
end

if isfield(opts,'tauD')
    if ~isinInterval(opts.tauD,0,2,false)
        error('If used, opts.tauD must be in (0,2).');
    end
end
    
if isfield(opts,'mxitr')
    if ~isinInterval(opts.mxitr,0,Inf,false) || opts.mxitr ~= floor(opts.mxitr)
        error('opts.mxitr should be a finite positive integer.');
    end
else
    opts.mxitr = 1000;
end

if isfield(opts,'eta')
    if ~isinInterval(opts.eta,1,Inf,false)
        error('opts.eta must be greater than one.');
    end
else
    opts.eta = 4;
end

if isfield(opts,'fullMu')
    if ~islogical(opts.fullMu)
        error('fullMu should be true or false.');
    end
else
    opts.fullMu = false;
end

if isfield(opts,'kappa')
    if ~isinInterval(opts.kappa,1,Inf,true)
        error('opts.kappa is a condition number and so should be >= 1.');
    end
end

if isfield(opts,'gamma')
    if ~isinInterval(opts.gamma,0,1,true)
        error('opts.gamma must be in [0,1].');
    end
else
    opts.gamma = 0.85;
end

if isfield(opts,'c')
    if ~isinInterval(opts.c,0,1,false)
        error('opts.c must be in (0,1).');
    end
else
    opts.c = 1E-3;
end

if isfield(opts,'beta')
    if ~isinInterval(opts.beta,0,1,false)
        error('opts.beta must be in (0,1).');
    end
else
    opts.beta = 0.5;
end

if isfield(opts,'scale')
    if ~islogical(opts.scale)
        error('scale should be true or false.');
    end
else
    opts.scale = true;
end

return