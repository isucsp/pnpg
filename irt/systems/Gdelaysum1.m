  function ob = Gdelaysum(varargin)
%|function ob = Gdelaysum([options])
%|
%| Construct Gdelaysum object, which performs weighted sums of a delayed signal.
%| This is useful for model-based reconstruction of THz images.
%| See Gdelaysum_test() at end for example usage.
%|
%| y[n;m] = sum_{k=0}^{Nx-1} x[k] h[n - d[k;m]], n=0,...,Ny-1, m=0,...,Nm-1
%| for 0 <= n - d[k] <= Nh - 1
%|
%| required
%|	'Ny'	[1]		output signal length
%|	'delay'	[Nx Nm]		delays
%|
%| options
%|	'h'	[Nh 1]		impulse response, default [1]
%|	'nthread' [1]		# threads, default 1
%|
%| out
%|	ob	[Ny*Nm Nx]	Fatrix object
%|
%| Copyright 2006-8-25, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), Gdelaysum_test, return, end
if nargin < 4, help(mfilename), error(mfilename), end

% defaults
arg.Ny = [];
arg.delay = [];
arg.h = 1;
arg.nthread = 1;

% options
arg = vararg_pair(arg, varargin);

if isempty(arg.delay), fail 'delay required', end
if isempty(arg.Ny), fail 'Ny required', end

arg.h = single(arg.h);
arg.delay = single(arg.delay);
arg.nthread = int32(arg.nthread);
arg.Ny = int32(arg.Ny);

[arg.Nx arg.Nm] = size(arg.delay);
arg.Nh = length(arg.h);

arg.odim = [arg.Ny*arg.Nm];
arg.dim = [arg.Ny * arg.Nm, arg.Nx];

%
% build Fatrix object
%
ob = Fatrix(arg.dim, arg, ...
	'forw', @Gdelaysum_forw, ...
	'back', @Gdelaysum_back, ...
	'caller', mfilename);


%
% Gdelaysum_forw(): y = G * x
% in
%	x	[Nx L]
% out
%	y	[Ny*Nm L]
%
function y = Gdelaysum_forw(arg, x)

LL = size(x, 2);
y = zeros(arg.Ny*arg.Nm, LL, 'single');
for ll=1:LL
	tmp = single(x(:,ll));
%keyboard
	tmp = delaysum1_mex('delaysum1,forw', arg.h, arg.delay, ...
		arg.nthread, tmp, arg.Ny);
	y(:,ll) = tmp(:);
end


%
% Gdelaysum_back(): x = G' * y
% in
%	y	[Ny*Nm L]
% out
%	x	[Nx L]
%
function x = Gdelaysum_back(arg, y)

LL = size(y, 2);
x = zeros(arg.Nx, LL, 'single');
for ll=1:LL
	tmp = single(y(:,ll));
	tmp = delaysum1_mex('delaysum1,back', arg.h, arg.delay, ...
		arg.nthread, tmp);
	x(:,ll) = tmp;
end


%
% Gdelaysum_test
%
function Gdelaysum_test

delay = [10 20 40; 50 51 52];
h = [4 3 2 1]';
Ny = 80;
A = Gdelaysum1('Ny', Ny, 'h', h, 'delay', delay);

x = [1 0]';
y1 = A * x;
y1 = reshapee(y1, Ny, []);
if im
	plot(y1, '-o')
end

mask = true(A.arg.Nx,1);
Fatrix_test_basic(A, mask)

test_adjoint(A);
