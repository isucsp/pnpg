% pwls_example.m
% example of how to use PWLS algorithms for PET reconstruction, such as PWLS-PCG
% Copyright 2003-11-30, Jeff Fessler, University of Michigan

% generate data
if ~isvar('yi'), printm 'setup'
	em_wls_test_setup
%	wi = ones(size(wi)); warning 'uniform wi' % to test circulant precon
	W = diag_sp(wi(:));
prompt
end

% regularization
if ~isvar('R'), printm 'R'
	f.l2b = 9;
	f.delta = 1;

	if 1
		Rq = Reg1(ig.mask, 'beta', 2^f.l2b);
		psf = qpwls_psf(G, Rq.C, 1, ig.mask);
		im(7, psf, 'PSF'), cbar
	end, clear Rq

	kappa = sqrt( div0(G' * wi, G' * sg.ones) );
	im(8, kappa, 'kappa'), cbar

	R = Reg1(kappa, 'type_denom', 'matlab', ...
		'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b);

	clear xos xiot xcg
prompt
end

f.niter = 20;

if ~isvar('xinit')
	%xinit = ones(size(xtrue));		% uniform initial image
	%xinit = xfbp;				% FBP initial image
	if exist('imfilter') == 2
		xinit = imfilter(xfbp, double(psf)); % smoothed FBP
	else
		xinit = convn(xfbp, double(psf), 'same'); % smoothed FBP
	end
	im(9, xinit, 'init'), cbar
prompt
end


if ~isvar('Gb'), printm 'do Gb'
	f.nsubset = 40;
	Gb = Gblock(G, f.nsubset);
prompt
end

% OS-SPS iterations (unconstrained)
if ~isvar('xos'), printm 'do os'
	xlim = [-inf inf]; % unconstrained
	[xos tim.os] = pwls_sps_os(xinit(ig.mask), yi, wi, Gb, R, ...
			1+f.niter, xlim, [], [], 1);
	xos = ig.embed(xos);
	im clf, im(xos, 'xos'), cbar
prompt
end

% PWLS-IOT iterations (unconstrained)
if 0 | ~isvar('xiot'), printm 'do iot'
	[xiot, tim.iot] = pl_iot(xinit(ig.mask), Gb, {yi, wi}, R, ...
			'dercurv', 'wls', ...
			'riter', 1, ...
			'os', 5, ... % f.niter, ...
			'niter', f.niter, 'isave', 'all', ...
			'pixmin', xlim(1), 'pixmax', xlim(2), ...
			'chat', 0);
	xiot = ig.embed(xiot);
	im clf, im(xiot, 'xiot'), cbar
	minmax(xiot-xos)
prompt
end

% CG iterations
if ~isvar('xcg'), printm 'xcg'
	[xcg, tim.cg] = pwls_pcg1(xinit(ig.mask), G, W, yi(:), R, ...
			'niter', f.niter, 'isave', 'all');
	xcg = ig.embed(xcg);
	im clf, im(xcg, 'CG'), cbar
prompt
end

% compare cost
if 1, printm 'cost plots'
	cost.cg		= pwls_cost(xcg,	G, W, yi(:), R, ig.mask);
	cost.os		= pwls_cost(xos,	G, W, yi(:), R, ig.mask);
	cost.iot	= pwls_cost(xiot,	G, W, yi(:), R, ig.mask);
	ii = 0:f.niter;
	if im
		clf, subplot(211)
		plot(ii, cost.os, 'y-o', ii, cost.cg, 'g-x', ii, cost.iot, 'c-+')
		xlabel 'iteration', ylabel 'cost'
		legend('OS-SPS', 'CG', 'IOT')

		subplot(212)
		plot(tim.os, cost.os, 'y-o', tim.cg(:,3), cost.cg, 'g-x', ...
			tim.iot, cost.iot, 'c-+')
		xlabel 'time', ylabel 'cost'
		legend('OS-SPS', 'CG', 'IOT')
	prompt
	end
%	minmax(diff(tim.os))
%	minmax(diff(tim.cg))
end

	pr nrms(col(xos(:,:,end)), xtrue(:))
	pr nrms(col(xcg(:,:,end)), xtrue(:))
	pr nrms(col(xiot(:,:,end)), xtrue(:))

if 1, printm 'images'
	clim = [0 6];
	im clf, im([xtrue, xfbp; xos(:,:,end), xiot(:,:,end); ...
		xcg(:,:,end), xinit], clim), cbar
	im plc 2 3
	im(1, xtrue, clim)
	im(2, xfbp, clim, 'FBP')
	im(3, xinit(:,:,end), clim, 'Init')
	im(4, xos(:,:,end), clim, 'OS')
	im(5, xiot(:,:,end), clim, 'IOT')
	im(6, xcg(:,:,end), clim, 'CG')
end
