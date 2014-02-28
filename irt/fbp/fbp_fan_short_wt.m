  function [wt scale180] = fbp_fan_short_wt(sg, varargin)
%|function [wt scale180] = fbp_fan_short_wt(sg, [options])
%|
%| Sinogram weighting for fan-beam short scan
%|
%| in
%|	sg	strum		sino_geom (sg.orbit_start is ignored)
%|
%| option
%|	'type'	'parker'	from parker:82:oss (Med Phys 1982)
%|
%| out
%|	wt	[nb na]		parker weights.  "excess" views have wt=0
%|	scale180		scale factor needed for back-projector	
%|
%| Usually you must multiply short fan-beam sinogram by
%| *both* the Parker weights "wt" and the scale factor
%| to get the proper results with my back projectors because
%| those back-projectors scale by the orbit and number of view angles.
%|
%| Copyright 2009-12-10, Jeff Fessler and Janghwan Cho, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(sg, 'test')
	[wt scale180] = fbp_fan_short_wt_test;
	if ~nargout, clear wt scale180, end
return
end

arg.type = 'parker';
arg = vararg_pair(arg, varargin);

switch arg.type
case 'parker'
	[wt scale180] = fbp_fan_short_wt_parker(sg);
otherwise
	fail('unknown type %s', arg.type)
end


% fbp_fan_short_wt_parker()
function [wt scale180] = fbp_fan_short_wt_parker(sg)

if sg.orbit < sg.orbit_short - eps
	warn('orbit %g is less than a short scan %g', sg.orbit, sg.orbit_short)
end

if sg.orbit > sg.orbit_short + eps
	warn('orbit %g exeeds short scan %g by %g views', ...
		sg.orbit, sg.orbit_short, ...
		(sg.orbit -  sg.orbit_short) / (sg.ad(2) - sg.ad(1)));
end

nb = sg.nb;
na = sg.na;
bet = sg.ar - sg.ar(1); % trick: force 0 start, so this ignores orbit_start!
gam = sg.gamma;
[gg bb] = ndgrid(gam, bet);
gam_max = sg.gamma_max; % half of fan angle

fun = @(x) sin(pi/2 * x).^2; % smooth out [0,1] ramp
% todo: could use integrals of this function over the
% tiny angular range of each projection view so that
% the sum over beta of these functions is a constant.
% or use a bspline (quadratic?) that would be easier to integrate?

wt = zeros(nb,na); % any extra views will get 0 weight
ii = bb < 2 * (gam_max - gg); % 0 <= bb not needed
tmp = bb(ii) ./ (2 * (gam_max - gg(ii)));
wt(ii) = fun(tmp);

ii = 2 * (gam_max - gg) <= bb & bb < pi - 2 * gg;
wt(ii) = 1;

ii = pi - 2 * gg < bb & bb <= pi + 2 * gam_max;
tmp = (pi + 2*gam_max - bb(ii)) ./ (2 * (gam_max + gg(ii)));
wt(ii) = fun(tmp);

scale180 = sg.orbit / 180; % scale factor due to orbit and na in backprojector


% fbp_fan_short_wt_test
function [wt scale180] = fbp_fan_short_wt_test

% generate object, image geometry
ig = image_geom('nx', 128, 'dx', 4);
ell = [	0 0 50 40 0 1;
	-18 0 12 12 0 1;
	18 -0 12 12 0 1] * 4;
xtrue = ellipse_im(ig, ell, 'oversample', 4);

% sinogram geometry for both short and full scan
orbit_start = 97; % stress test
sg_short = sino_geom('ge1', 'orbit', 'short', 'down', 4, ...
	'orbit_start', orbit_start);
if 1 % examine effect of extra views
	dd = sg_short.ad(2) - sg_short.ad(1);
	extra = 10;
	sg_short.na = sg_short.na + extra;
	sg_short.orbit = sg_short.orbit + extra * dd;
end

sg_360 = sino_geom('ge1', 'down', sg_short.down, ...
	'orbit_start', sg_short.orbit_start);

% sinogram of the object for each scan
sino_short = ellipse_sino(sg_short, ell, 'oversample', 2);
sino_360 = ellipse_sino(sg_360, ell, 'oversample', 2);

% apply parker weighting
[wt scale180] = fbp_fan_short_wt(sg_short); % parker weighting
sino_parker = sino_short .* wt;

% FBP reconstructed images
fbp_geom_short = fbp2(sg_short, ig);
fbp_geom_360 = fbp2(sg_360, ig);

fbp_w_short = fbp2(sino_short, fbp_geom_short);
fbp_w_parker = scale180 * fbp2(sino_parker, fbp_geom_short);
fbp_w_360 = fbp2(sino_360, fbp_geom_360);

% plot
im plc 4 3
im subplot [1 2 3]
im(rad2deg(sg_short.gamma), sg_short.ad, wt), cbar
xlabel 'gamma [degrees]'
ylabel 'beta [degrees]'
title 'Parker weighting'
clim = [0 8];
im(4, fbp_w_360, clim, 'FBP: full scan'), cbar
im(5, fbp_w_short, clim, 'FBP: short scan w/o parker weighting'), cbar
im(6, fbp_w_parker, clim, 'FBP, short scan w/ parker weighting'), cbar
im(7, xtrue, clim, 'True'), cbar
im subplot [8 9]
plot([fbp_w_360(:,end/2) fbp_w_short(:,end/2) fbp_w_parker(:,end/2)])
title 'middle slice, y = end/2'
legend('full', 'short w/o parker', 'short w/ parker')

im(10, sino_360)
im(11, sino_short)
im(12, sino_parker)

%clf, im_toggle(fbp_w_360 .* ig.circ, fbp_w_parker .* ig.circ)
%clf, plot(scale180 * mean(wt'), '.') % within 0.005 of 1

%yaxis_pi '0 p'
%plot(diff(wt,1))
%savefig fig_tomo_fan_short_wt
