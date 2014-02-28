  function ob = Gdft(varargin)
%|function ob = Gdft([args])
%| Construct Gdft object that computes samples of the DFT of a signal
%| with dimensions [(Nd)].  This is useful for "under-sampled" MRI.
%|
%| options (at least one of first two must be provided):
%|	'mask'	logical [(Nd)]	image-domain mask, usually: true(nx,ny)
%|	'samp'	logical [(Nd)]	which frequency-domain samples to return
%|				(default: all samples)
%|	'ifftshift'	0|1	apply ifftshift to image before DFT?
%|	'fftshift'	0|1	apply fftshift to spectrum after DFT?
%|				(both default to false)
%|
%| out
%|	ob	[M np]	Fatrix object, where M = sum(samp(:)), np = sum(mask(:))
%|
%| See Gdft_test.m for example usage.
%|
%| Basically, you create a system matrix object by calling:
%|	A = Gdft( ... )
%| and then you can use it thereafter by typing commands like
%|	y = A * x;
%| which will auto-magically evaluate the DFT samples.
%| This is useful for iterative image reconstruction in MRI.
%|
%| Besides simple utilities like display, there are the following
%| capabilities of this object:
%|	y = A * x		forward operation
%|	x = A' * y		adjoint operation
%|	
%| Copyright 2003-6-1, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(varargin{1}, 'test'), Gdft_test, return, end

arg.mask = [];
arg.samp = [];
arg.ifftshift = false; % apply ifftshift to image before DFT?
arg.fftshift = false; % apply fftshift to spectrum after DFT?
arg = vararg_pair(arg, varargin);

if isempty(arg.mask) && isempty(arg.samp), error 'must give mask or samp', end

if isempty(arg.mask)
	arg.mask = true(size(arg.samp));
elseif isempty(arg.samp)
	arg.samp = true(size(arg.mask));
else
	if ~isequal(size(arg.mask), size(arg.samp)), fail 'mismatch', end
end

if ~islogical(arg.mask), error 'mask must be logical', end
if ~islogical(arg.samp), error 'samp must be logical', end

arg.Nd = size(arg.mask);
arg.ndim = length(arg.Nd);
if arg.Nd(end) == 1 % 1D case
	arg.ndim = arg.ndim - 1;
end
arg.M = sum(arg.samp(:));
arg.np = sum(arg.mask(:));
arg.dim = [arg.M arg.np];

% build Fatrix object
ob = Fatrix(arg.dim, arg, ...
	'forw', @Gdft_forw, 'back', @Gdft_back, ...
	'gram', @Gdft_gram, 'caller', mfilename);


%
% Gdft_forw(): y = G * x
% in
%	x	[np L] or [(Nd) L]
% out
%	y	[M L]
%
function y = Gdft_forw(arg, x)

if size(x,1) == arg.np		% [np (L)]
	x = embed(x, arg.mask);	% [(Nd) (L)]
end

if isequal(size(x), size(arg.mask)) % [(Nd)]
	y = fftn_shifted(x, arg.ifftshift, arg.fftshift); % [(Nd)]
	y = y(arg.samp); % [M 1]
else
	xdim = size(x);
	NN = prod(arg.Nd);
	LL = xdim(end);

	if arg.ndim == 1
%		y = fft(x, [], 1);
		y = fft1_shifted(x, arg.ifftshift, arg.fftshift); % [Nd (L)]
		y = y(arg.samp,:); % [M (L)]

	else
		y = zeros([NN LL]);
		for ll=1:LL
			tmp = stackpick(x, ll);
			y(:,ll) = col(fftn_shifted(tmp, ...
				arg.ifftshift, arg.fftshift));
		end
		y = y(arg.samp,:); % [M L]
	end
end


%
% Gdft_back(): x = G' * y
% in
%	y	[M L]
% out
%	x	[np L]
%
function x = Gdft_back(arg, y)

y = embed(y, arg.samp); % [(Nd) L]
ydim = size(y);
NN = prod(arg.Nd);

if isequal(size(y), size(arg.samp)) % [(Nd)]
	x = ifftn_shifted(y, arg.ifftshift, arg.fftshift); % [(Nd)] adjoint
	x = col(x); % [*Nd]
else
	if arg.ndim == 1
		x = ifft1_shifted(y, arg.ifftshift, arg.fftshift); % [Nd L]
	else
		LL = ydim(end);
		x = zeros([NN LL]);
		for ll=1:LL
			tmp = stackpick(y,ll);
			x(:,ll) = col(ifftn_shifted(tmp, ...
				arg.ifftshift, arg.fftshift));
		end
	end
end

% note the "NN" factor that is needed to make it an adjoint operation:
x = NN * x(arg.mask(:),:); % [np *L]

%x = reshape(x, [Ns Ld]); % [np (L)] % not needed if L can be scalar only!


%
% fft1_shifted()
%
function y = fft1_shifted(x, is_ifftshift, is_fftshift)
if is_ifftshift
	x = ifftshift(x, 1);
end
y = fft(x, [], 1);
if is_fftshift
	y = fftshift(y, 1);
end


%
% ifft1_shifted()
%
function x = ifft1_shifted(y, is_ifftshift, is_fftshift)
if is_fftshift
	y = fftshift(y, 1);
end
x = ifft(y, [], 1);
if is_ifftshift
	x = ifftshift(x, 1);
end


%
% fftn_shifted()
%
function y = fftn_shifted(x, is_ifftshift, is_fftshift)
if is_ifftshift
	x = ifftshift(x);
end
y = fftn(x);
if is_fftshift
	y = fftshift(y);
end


%
% ifftn_shifted()
%
function x = ifftn_shifted(y, is_ifftshift, is_fftshift)
if is_fftshift
	if any(rem(size(y), 2))
		error 'odd size not done with shift'
	end
	y = fftshift(y);
end
x = ifftn(y);
if is_ifftshift
	if any(rem(size(x), 2))
		error 'odd size not done with shift'
	end
	x = ifftshift(x);
end


%
% Gdft_gram()
%
function [T, reuse] = Gdft_gram(A, W, reuse)
if isempty(W)
	T = A' * A;
else
	T = A' * W * A;
end


% Gdft_test
function Gdft_test

Nd = [10 6];
rand('state', 1)
samp = rand(Nd) > 0.3;
mask = true(Nd); mask(1) = false; % stress

for ii = [1 0]
for jj = [1 0]
	A = Gdft('mask', mask, 'samp', samp, 'ifftshift', ii, 'fftshift', jj);
	Fatrix_test_basic(A, mask, 'complex', 1)
	test_adjoint(A, 'complex', 1, 'tole', 2e-16);
end
end

if 1 % test vs fft
	x = rand(Nd);
	x = x .* mask;
	y = fftn(x);
	y = y(samp(:));
	jf_equal(y, A * x)

	x = ifftn(embed(y, samp)) * prod(Nd);
	x = x .* mask;
	jf_equal(x, embed(A' * y, mask))
end

if 1 % test gram
	wi = col(1:sum(samp(:)));
	W = diag_sp(wi);
	T = build_gram(A, W);
	y1 = T * x;
	y2 = (A' * (wi .* (A * x)));
	jf_equal(y1, y2)
end
