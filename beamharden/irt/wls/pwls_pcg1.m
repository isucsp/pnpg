  function [xs, info] = pwls_pcg1(x, A, W, yi, R, varargin)
%|function [xs, info] = pwls_pcg1(x, A, W, yi, R, [options])
%|
%| penalized weighted least squares (PWLS)
%| with convex non-quadratic regularization,
%| minimized via preconditioned conjugate gradient algorithm.
%| cost(x) = (y-Ax)'W(y-Ax)/2 + R(x)
%| See pwls_example.m for usage.
%|
%| in
%|	x	[np 1]		initial estimate
%|	A	[nd np]		system matrix
%|	W	[nd nd]		data weighting matrix, usually diag_sp(wi)
%|	yi	[nd 1]		noisy data
%|	R			penalty object (see Robject.m)
%|
%| options
%|	niter			# total iterations
%|	isave	[]		list of iterations to archive (default: 'last')
%|	userfun			user defined function handle (see default below)
%|	precon	[np np]		preconditioner (default: 1, i.e., none)
%|	stepper			method for step-size line search
%|				default: {'qs', 3}
%|	stop_threshold		stop iterations if norm(xnew-xold)/norm(xnew)
%|				is less than this unitless value.  default: 0
%|	stop_norm_type		use norm(.,type) for stop rule
%|				choices: 1 | 2 (default) | inf
%|	chat	0|1		verbosity (default 0)
%|
%| out
%|	xs	[np niter]	estimates each iteration
%|	info	[niter 3]	gamma, step, time
%|
%| Copyright 1996-7, Jeff Fessler, University of Michigan

if nargin == 1 && streq(x, 'test'), pwls_pcg1_test, return, end
if nargin < 5, help(mfilename), error args, end

% defaults
arg.precon = 1;
arg.niter = 1;
arg.isave = [];
arg.stepper = {'qs', 3};	% quad surr with this # of subiterations
arg.userfun = @userfun_default;
arg.stop_threshold = 0;
arg.stop_norm_type = 2;
arg.chat = false;
arg.key = 1;

arg = vararg_pair(arg, varargin);

arg.isave = iter_saver(arg.isave, arg.niter);
mynorm = @(x) norm(x, arg.stop_norm_type);

cpu etic

x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(niter,?); % trick: do not initialize since size may change

%
% initialize projections
%
ticker(mfilename, 1, arg.niter)
Ax = A * x;

oldinprod = 0;

%
% iterate
%
for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	%
	% (negative) gradient
	%
	ngrad = A' * (W * (yi-Ax));
	pgrad = R.cgrad(R, x);
	ngrad = ngrad - pgrad;

	%
	% preconditioned gradient
	%
	pregrad = arg.precon * ngrad;

	%
	% direction
	%
	newinprod = ngrad' * pregrad;
	newinprod = reale(newinprod, 'warn', 'inprod');
	if iter == 1
		ddir = pregrad;
		gamma = 0;
	else
		if oldinprod == 0
			warning 'inprod=0.  going nowhere!'
			gamma = 0;
		else
			gamma = newinprod / oldinprod;	% Fletcher-Reeves
%			gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = ngrad;
	oldinprod = newinprod;

	% check if descent direction
	if real(ddir' * ngrad) < 0
		warning 'wrong direction'
		if arg.key, keyboard, end
	end

	%
	% step size in search direction
	%
	Adir = A * ddir;
%	Cdir = R.C * ddir; % this is too big for 3D CT problems

	%
	% one step based on quadratic surrogate for penalty
	%
	if streq(arg.stepper{1}, 'qs1')
%		pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); % avoid Cdir
		pdenom = (abs(ddir).^2)' * R.denom(R, x);
		denom = Adir'*(W*Adir) + pdenom;
		if denom == 0
			warning 'found exact solution???  step=0 now!?'
			step = 0;
		else
			step = real((ddir' * ngrad) / denom);
		end

	%
	% Iteratively minimize the 1D line-search function over step:
	%	1/2 || y - A (x + step*ddir) ||_W^2 + R(x + step*ddir)
	% This is a real-valued function of the real-valued step parameter.
	%
	elseif streq(arg.stepper{1}, 'qs')
		nsub = arg.stepper{2};
		dAWAd = Adir' * (W * Adir);
		dAWAd = reale(dAWAd); % 2008-10-16
		dAWr = Adir' * (W * (yi-Ax));
		dAWr = real(dAWr); % 2008-10-16
		step = 0;
		for is=1:nsub
%			pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); % avoid Cdir
			pdenom = (abs(ddir).^2)' * R.denom(R, x + step * ddir);
			denom = dAWAd + pdenom;
			if denom == 0 || isinf(denom)
				'0 or inf denom?'
				if arg.key, keyboard, end
				error bad
			end
			pgrad = R.cgrad(R, x + step * ddir);
			pdot = real(ddir' * pgrad); % 2008-10-15
			step = step - (-dAWr + step * dAWAd + pdot) / denom;
%			2008-10-16: removed below because made real above
%			step = real(step); % real step size seems logical
		end

	else
		error 'bad stepper'
	end

	if step < 0
		warning 'downhill?'
		if arg.key, keyboard, end
	end

	%
	% update
	%
	Ax = Ax + step * Adir;
%	Cx = Cx + step * Cdir;
	x = x + step * ddir;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
        end
        info(1+iter,:) = feval(arg.userfun);

	% check norm(xnew-xold) / norm(xnew) vs threshold
	if arg.stop_threshold && ...
		mynorm(step * ddir) / mynorm(x) < arg.stop_threshold
		if arg.chat
			printm('stop at iteration %d', iter)
		end
		if isequal(arg.isave, arg.niter) % saving last iterate only?
			xs(:,1) = x; % be sure to save this 'final' iterate
		else % saving many iterates?
			xs(:, arg.isave > iter) = []; % clear out unused
		end
		return
	end
end

% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default
gamma = evalin('caller', 'gamma');
step = evalin('caller', 'step');
out = [gamma step cpu('etoc')];



function pwls_pcg1_test
mask = true([11 12]); mask(1) = false;
A = 2 * Gdft('mask', mask); % orthogonal to permit analytical solution
xtrue = zeros(size(mask), 'single');
xtrue(round(end/2), end/2) = 1;
y = A * xtrue;
beta = 2^6;
R = Reg1(mask, 'beta', beta);
xhat = full(A' * A + R.C' * R.C) \ (A' * y);
xhat = embed(xhat, mask);
im(xhat)

xinit = 0 * mask;
xpcg = pwls_pcg1(xinit(mask(:)), A, 1, y, R, 'niter', 100, ...
	'chat', 1, 'stop_threshold', 1e-3, 'stop_norm_type', 1);
xpcg = embed(xpcg, mask);

im clf, im pl 2 2
im(1, xtrue)
im(2, xhat)
im(3, xpcg)
im(4, xpcg - xhat)
