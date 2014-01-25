% interp_from_sparse.m
% An example of using PWLS-SPS algorithm to interpolate from sparsely
% sampled data.
% Copyright 2010-1-7, Jeff Fessler, University of Michigan

%
% generate data
%
if ~isvar('yi')
	ig = image_geom('nx', 110, 'ny', 128, 'dx', 1);
	xtrue = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2);
	im(221, xtrue, 'xtrue (high res)')

	% generate sparse "down sampling" sampling pattern
	samp = false(ig.nx, ig.ny);
	down = 2;
	samp(1:3:end,1:3:end) = true;

	yi = samp .* xtrue;
	im(222, yi, 'yi (sampled)')
prompt
end

%
% elementary 0-order hold zooming
% and cubic interpolation
%
if 0
if ~isvar('xc')
	x0 = down^2 * ig.shape(G' * yi(:));
	im(223, x0, 'x0, 0-order hold')

	xc = interp2(yi', 0.5+[0:(n.b*down-1)]'/down, 0.5+[0:(n.a*down-1)]/down)';
	xc(isnan(xc)) = 0;
	im(224, xc, 'xc, cubic interpolation')
prompt
end
end

% regularization
if ~isvar('Rn')
	Rq = Robject(ig.mask, 'type_denom', 'matlab', ...
		'potential', 'quad', 'beta', 2^(-4));
	Rn = Reg1(ig.mask, 'type_denom', 'matlab', ...
		'pot_arg', {'huber', 0.3}, 'beta', 2^(-7));
prompt
end

%
% SPS with quadratic penalty
%
W = diag_sp(samp(:));
xinit = yi;
if ~isvar('xq')
	f.niter = 40;
	A = diag_sp(ones(ig.nx*ig.ny,1));
	xq = pwls_sps_os(yi(:), yi(:), W, A, Rq, f.niter, inf, [], [], 1, 0);
	xq = ig.embed(xq);

	im clf, im(xq, 'QPWLS-SPS iterations')
prompt
end

return

%
% SPS with nonquadratic edge-preserving penalty
%
if ~isvar('xs'), printm 'xs'
	f.niter = 50-1;
	xs = pwls_sps_os(col(xq(:,:,end)), yi, [], G, Rn, ...
		f.niter, 10, [], [], 1, 0);
	xs = ig.embed(xs);

	im clf, im(xs, 'PWLS-SPS iterations')
prompt
end

if im
	rms = inline('sqrt(mean(abs(x).^2))', 'x');
	r0 = rms(col(xtrue-x0));
	rc = rms(col(xtrue-xc));
	rq = rms(col(xtrue-xq(:,:,end)));
	rs = rms(col(xtrue-xs(:,:,end)));

	im pl 2 3
	clim = [0 6];
	im(1, xtrue, clim, 'truth'), cbar
	im(2, yi, clim, 'yi'), cbar
	axis([1 ig.nx 1 ig.ny])
	im(3, x0, clim, 'x0, zero-order'), cbar
	xlabel(num2str(r0))
	im(4, xc, clim, 'xc, cubic'), cbar
	xlabel(num2str(rc))
	im(5, xq(:,:,end), clim, 'SPS quadratic'), cbar
	xlabel(num2str(rq))
	im(6, xs(:,:,end), clim, 'SPS edge preserving'), cbar
	xlabel(num2str(rs))
end
