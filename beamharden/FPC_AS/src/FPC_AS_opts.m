function opts = FPC_AS_opts(opts,n,m)
% Options for FPC_ActiveSet 
%
%--------------------------------------------------------------------------
% DESCRIPTION
%--------------------------------------------------------------------------
%
% opts = fpc_active_opts(opts)
%
% If opts is empty upon input, opts will be returned containing the default options  
%
% Alternatively, if opts is passed with some fields already defined, those
% fields will be checked for errors, and the remaining fields will be added
% and initialized to their default values.
%
% (This file should be cleaned. Some of fields are not used at all.)
%
% Copyright (c) 2008.  Zaiwen Wen, Wotao Yin
%
%--------------------------------------------------------------------------
%%
%
%
if isempty(opts); opts.memo =[]; end
% intial value
opts = AddLowerBoundedNumberOption(opts, 'x0', [], -inf, inf, 'initial solution');
opts = AddLowerBoundedNumberOption(opts, 'init', 2, 0, 10, 'methods of initialization, integer');
opts = AddLowerBoundedNumberOption(opts, 'xs', [], -inf, inf, 'exact solution whose purpose is only for comparison');
opts = AddLowerBoundedNumberOption(opts, 'tol_eig', 1e-4, 0, 1, 'tolerance for eigs');
opts = AddLowerBoundedNumberOption(opts, 'scale_A', 0, 0, 1, 'scale the matrix A so that max of eigs(A*AT) equals 1, integer');
opts = AddLowerBoundedNumberOption(opts, 'eigs_mxitr', 20, 1, 100, 'max number of iterations for eigs(A*AT), integer');


% estimation of the solution
opts = AddLowerBoundedNumberOption(opts, 'eps', 1e-16, 0, 1,'machine accuarcy');
opts = AddLowerBoundedNumberOption(opts, 'zero', 1e-10, 0, 1e10,'minimal magnitude of x');
opts = AddLowerBoundedNumberOption(opts, 'dynamic_zero', 1, 0, 1,'set the thresholding level dynamically or not, integer');
opts = AddLowerBoundedNumberOption(opts, 'minK', floor(m/2), 1, n,'etimation of the minimal number of the nonzero components of the exact solution, integer');
opts = AddLowerBoundedNumberOption(opts, 'maxK', m, 1, n,'etimation of the maximal number of the nonzero components of the exact solution, integer');
opts = AddLowerBoundedNumberOption(opts, 'hard_truncate', 1, 0, 1,'do hard truncation before doing sub-opt, integer');


% options for shrinkage
opts = AddLowerBoundedNumberOption(opts, 'tauD', min(1.999,-1.665*m/n + 2.665), 0, 1e5, 'a parameter for shrinkage');
opts = AddLowerBoundedNumberOption(opts, 'tau_min', 1e-4, 0, 1e5, 'minimal tau');
opts = AddLowerBoundedNumberOption(opts, 'tau_max', 1e3, 0, 1e5, 'minimal tau');

% options for stopping rules
opts = AddLowerBoundedNumberOption(opts, 'mxitr', 1000, 1, 1e5, 'max number of iterations, integer');
%opts = AddLowerBoundedNumberOption(opts, 'xtol', 1e-8, 0, 1, 'Tolerance on norm(x - xp)/norm(xp)');
opts = AddLowerBoundedNumberOption(opts, 'gtol', 1e-6, 0, 1,'Tolerance on norm of sub-gradient');
opts = AddLowerBoundedNumberOption(opts, 'gtol_scale_x', 1e-12, 0, 1,'Tolerance on norm of sub-gradient');
opts = AddLowerBoundedNumberOption(opts, 'f_rel_tol', 1e-20, 0, 1,'Tolerance on the relative changes of function value');
opts = AddLowerBoundedNumberOption(opts, 'f_value_tol', 0, 0, inf,'Tolerance on the optimal objective value');

% options for line search
opts = AddLowerBoundedNumberOption(opts, 'ls_mxitr', 5, 2, 100,'max number of iterations of of line search subroutine, integer');
opts = AddLowerBoundedNumberOption(opts, 'gamma', 0.85, 0, 1,'a parameter for the nonmonotone line search');
opts = AddLowerBoundedNumberOption(opts, 'c', 1e-3, 0, 1,'a parameter for Armijo condition');
opts = AddLowerBoundedNumberOption(opts, 'beta', 0.5, 0, 1,'a parameter for decreasing the step size in the nonmonotone line search');

% options for continuation
opts = AddLowerBoundedNumberOption(opts, 'eta',  0.1, 0, 1,'a parameter for decreasing mu');
% opts = AddLowerBoundedNumberOption(opts, 'eta_rho',  0.5, 0, 1,'a parameter for decreasing eta');
% opts = AddLowerBoundedNumberOption(opts, 'eta_min',  0.01, 0, 1,'minimal eta');
% opts = AddLowerBoundedNumberOption(opts, 'eta_max',  0.5, 0, 1,'maximal eta');

% options for sub-optimization
opts = AddLowerBoundedNumberOption(opts, 'sub_mxitr', 50, 1, 1e5, 'max number of iterations for doing sub-optimization, integer');
opts = AddLowerBoundedNumberOption(opts, 'lbfgs_m', 5, 1, 100, 'storage number of L-BFGS, integer');
opts = AddStringOption(opts, 'ls_meth', 'Hybridls', {'nmls','Hybridls'}, 'line search methods, string');
opts = AddStringOption(opts, 'sub_opt_meth', 'pcg', {'lbfgs','lbfgsb','pcg'}, 'sub-optimization methods, string');

opts = AddLowerBoundedNumberOption(opts, 'kappa_g_d', 10, 1, 1e3,'tolerance for checking whether do sub-optimization or not');
opts = AddLowerBoundedNumberOption(opts, 'kappa_rho', 10, 1, 1e3,'a parameter for increasing kappa_g_d');

opts = AddLowerBoundedNumberOption(opts, 'tol_start_sub', 1e-6, 0, 1,'Tolerance of starting sub-optimization');

opts = AddLowerBoundedNumberOption(opts, 'min_itr_shrink', 3, 1, 1e3, 'min number of itertations between two sub-optimization, integer');
opts = AddLowerBoundedNumberOption(opts, 'max_itr_shrink', 20, 1, 1e3,'max number of itertations between two sub-optimization, integer');


% options for debugging
%opts = AddLowerBoundedNumberOption(opts, 'plot', 0, 0, 1,'plot solution');
opts = AddLowerBoundedNumberOption(opts, 'record', 0, -100, 100,'print information, -1=quiet, 0=some output, 1=more output. integer');
opts = AddLowerBoundedNumberOption(opts, 'PrintOptions', 0, 0, 1,'print options, integer');

end

