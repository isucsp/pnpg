% cone_beam_ct_example.m
% Illustrate cone-beam X-ray CT image reconstruction via FDK and iterative.
% This illustration uses a tiny system size cone-beam projection can be slow.
% Copyright 2005-1-21, Jeff Fessler, University of Michigan

% First run the FDK example.  It generates the true image xtrue
% and noiseless projection views "proj" and noisy data "yi"
% and generates (noisy) FDK recon "xfdk" for comparison / initialization.
if ~isvar('xfdk')
	bi = 1e6; % 1M photons / ray
	ri = 0; % no scatter etc. for now
	dfs = inf; % flat!
	feldkamp_example
prompt
end


% system matrix
if ~isvar('A'), printm 'A'
	f.nz_pad = 18; % add this many slices to each side.  todo: how many?
	ig_pad = ig.expand_nz(f.nz_pad)
	A = Gcone(cg, ig_pad, 'type', 'sf1');
%	f.sys_type = aspire_pair(cg, ig_pad, 'system', '3l'); % old way
%	A = Gtomo3(f.sys_type, ig_pad.mask, ig.nx, ig.ny, ig_pad.nz, ...
%		'chat', 0, 'permute213', true, 'checkmask', im&0);
end


% block object for ordered-subsets iterations
if ~isvar('Ab'), printm 'Ab'
	f.nblock = 8;
	Ab = Gblock(A, f.nblock);
end


if 1 % padding functions to help with long object problem
	f.rep = @(x) cat(3, repmat(x(:,:,1), [1 1 f.nz_pad]), x, ...
		repmat(x(:,:,end), [1 1 f.nz_pad]));
	f.pad = @(x) cat(3, zeros(ig.nx, ig.ny, f.nz_pad), x, ...
		zeros(ig.nx, ig.ny, f.nz_pad));
end


% check discrete vs analytical projection (note long object problem)
if 0, printm 'proj check'
	cpu etic
	pp = Ab * f.rep(xtrue);
	cpu etoc 'proj time'
	nrms(pp, proj)
	im clf, im_toggle(proj(:,:,1:12:end), pp(:,:,1:12:end), [0 4.4])
prompt
end


% regularization object
if ~isvar('R'), printm 'regularizer'
	f.l2b = 2^4.5;
	f.delta = 100/1000;
	R = Reg1(ig_pad.mask, 'type_denom', 'matlab', ...
                'pot_arg', {'hyper3', f.delta}, 'beta', 2^f.l2b);
	if 1 % check spatial resolution (away from edges)
		W = diag_sp(yi(:));
		psf = qpwls_psf(A, R, 1, ig_pad.mask, W, 'fwhmtype', 'profile');
	prompt
	end
end

	xinit = ig_pad.maskit(f.rep(xfdk));

if ~isvar('xpwls'), printm 'PWLS reconstruction'
	si = log(bi ./ yi); % log sinogram
	wi = yi;
	xpwls = pwls_pcg1(xinit, A, diag_sp(wi), si(:), R, ...
		'niter', 100, 'stop_threshold', 1e-3);
	xpwls = ig_pad.embed(xpwls);
	im(xpwls), cbar
prompt
end


% reshape data to be "2d arrays" for OS iterations (subset over last dim)
if ~isvar('os_data'), printm 'os_data'
	if isscalar(bi) && isscalar(ri)
		os_data = reshaper(yi, '2d');
		os_data = {os_data, ...
			bi * ones(size(os_data)), ri * ones(size(os_data))};
	else
		os_data = {reshaper(yi, '2d'), reshaper(bi, '2d'), ...
			reshaper(ri, '2d')}; % all data as 2d arrays
	end
end


% OS-SPS iterations for transmission penalized likelihood
if ~isvar('xpl'), printm 'start iterations'
	f.niter = 20;
	xs = tpl_os_sps(xinit, Ab, os_data{:}, R, 1+f.niter);
	xs = ig_pad.embed(xs);
	xpl = xs(:,:,:,end);
	im(xpl)
end


% finally, compare FDK vs iterative results
if 1
	iz_good = f.nz_pad + [1:ig.nz]; % original slices
	xpwls_good = xpwls(:,:,iz_good);
	xpl_good = xpl(:,:,iz_good);
	clim = [0 0.02];
	im plc 2 2
	im(1, xtrue, 'true', clim), cbar
	im(2, xfdk, 'FDK', clim), cbar
	im(4, xpwls_good, 'PWLS', clim), cbar
	im(3, xpl_good, 'PL', clim), cbar
	nrms(xtrue, xfdk)
	nrms(xtrue, xpl_good)
	nrms(xtrue, xpwls_good)
end
