  function [sino, pos, ang] = rect_sino(sg, rects, varargin)
%|function [sino, pos, ang] = rect_sino(sg, rects, [options])
%|
%| Create sinogram projection of one or more rectangles.
%| Works for both parallel-beam geometry and for fan-beam geometry,
%| i.e., 3rd generation x-ray ct systems, if the distances are specified.
%|
%| in:
%|	sg			sinogram geometry object from sino_geom()
%|	rects	[ne 6]		[centx centy widthx widthy angle_degrees amplitude]
%|
%| options:
%|	'oversample'		oversampling factor for emulating "strips"
%|	'xscale'		use -1 to flip in x (not recommended)
%|	'yscale'		use -1 to flip in y (not recommended)
%|
%| out:
%|	sino	[nb na]		sinogram
%|	pos	[nb 1]		radial sample positions (or [nb na] if mojette)
%|	ang	[na]		angular sample locations
%|
%| 2008-08-04, Yong Long, adapted from ellipse_sino()
%| Copyright 2003-10-22, Jeff Fessler, University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if streq(sg, 'test')
	rect_sino_test % no test on moj
return
end

% defaults
arg.xscale = 1;
arg.yscale = 1;
arg.oversample = 1;
arg = vararg_pair(arg, varargin);

if streq(sg.type, 'fan')
	[sino pos ang] = rect_sino_go(rects, [], [], ...
	sg.nb, sg.ds, sg.offset_s, ...
	sg.na, sg.orbit, sg.orbit_start, ...
	sg.dso, sg.dod, sg.dfs, ...
	sg.source_offset, ...
	arg.xscale, arg.yscale, arg.oversample, 0);

elseif streq(sg.type, 'par')
	[sino pos ang] = rect_sino_go(rects, [], [], ...
	sg.nb, sg.dr, sg.offset_r, ...
	sg.na, sg.orbit, sg.orbit_start, ...
	inf, 1, 0, ...
	sg.source_offset, ...
	arg.xscale, arg.yscale, arg.oversample, 0);

elseif streq(sg.type, 'moj')
	[sino pos ang] = rect_sino_go(rects, [], [], ...
	sg.nb, [], sg.offset_r, ...
	sg.na, sg.orbit, sg.orbit_start, ...
	inf, 1, 0, ...
	sg.source_offset, ...
	arg.xscale, arg.yscale, arg.oversample, sg.dx);

else
	error 'sino geom not done'
end


%
% rect_sino_go()
%
function [sino, pos, ang] = rect_sino_go(rects, ...
	pos, ang, ...
	nb, ds, offset_s, ...
	na, orbit, orbit_start, ...
	dso, dod, dfs, ...
	source_offset, ...
	xscale, yscale, ...
	oversample, mojette);

if size(rects, 2) ~= 6, error '6 parameters per rect', end

if isempty(ang)
	ang = deg2rad(orbit_start + [0:na-1]/na * orbit);
end

[pos pos2] = rect_sino_pos(pos(:), nb, ds, offset_s, ...
		oversample, mojette, ang);

sino = rect_sino_do(rects, pos2, ang(:)', ...
	xscale, yscale, ...
	dso, dod, dfs, source_offset);

if oversample
	sino = downsample2(sino, [oversample 1]);
end


%
% rect_sino_pos()
% determine usual and fine "radial" sample positions
%
function [pos_coarse, pos_fine] = ...
	rect_sino_pos(pos, nb, ds, offset_s, nover, mojette, ang)

if isempty(pos)
	if isempty(nb), warning 'nb required when no pos provided', end
	wb = (nb-1)/2 + offset_s;
else
	if ~isempty(nb), warning 'nb ignored when pos provided', end
end

if mojette ~= 0 % tricky mojette radial sampling
	% trick: ray_spacing aka ds comes from dx which is in mojette
	if ~isempty(ds), error 'ds must be empty for mojette case', end
	dt = abs(mojette) * max(abs(cos(ang)), abs(sin(ang)))'; % [1,na]

	if ~isempty(pos), error 'mojette requires empty "pos"', end

	na = length(ang);
	pos_coarse = ([0:nb-1]-wb)' * dt(:)'; % [nb na]
	if nover > 1
		tmp = [-(nover-1):2:(nover-1)]' / (2*nover);
		pos_fine = zeros(nover*nb,na);
		for ia=1:na
			pos_fine(:,ia) = col(outer_sum(tmp, pos_coarse(:,ia)'));
		end
	else
		pos_fine = pos_coarse;
	end

else % ordinary non-mojette sampling

	if isempty(pos)
		pos_coarse = ([0:nb-1]'-wb) * ds; % [nb 1]
	else
		pos_coarse = pos(:); % [nb 1]
		ds = pos(2) - pos(1);
		if any(abs(diff(pos) / ds - 1) > 1e-10)
			error 'uniform spacing required'
		end
	end

	if nover > 1
		% determine fine sampling positions
		pos_fine = [-(nover-1):2:(nover-1)]' / (2*nover);
		pos_fine = outer_sum(pos_fine, pos_coarse(:)'); % [nover nb]
		pos_fine = pos_fine(:); % [nb*nover 1]
	else
		pos_fine = pos_coarse;
	end
end


%
% rect_sino_do()
%
function sino = rect_sino_do(rects, pos, ang, ...
	xscale, yscale, dso, dod, dfs, source_offset)

nb = size(pos,1);
na = length(ang);

%
% effective radial and angular sample locations in parallel geometry!
%
if isinf(dso)
	if size(pos,2) > 1 % mojette
		rads = pos;
		angs = repmat(ang, [nb 1]);	% [nb na]
	else
		[rads angs] = ndgrid(pos, ang);	% [nb na]
	end
else
	if size(pos,2) > 1, error 'mojette fan not supported', end
	dis_src_det = dso + dod;

	if isinf(dfs)	% flat detector
		gam = atan(pos / dis_src_det); % "gamma"
	else			% arc detector
		dis_foc_det = dfs + dis_src_det;
		alf = pos / dis_foc_det;
		gam = atan2(	dis_foc_det * sin(alf), ...
				dis_foc_det * cos(alf) - dfs );	% gamma
	end

	rad = dso * sin(gam) + source_offset * cos(gam);	
	rads = rad(:,ones(1,na));
	angs = outer_sum(gam, ang); % gamma + beta
end
clear alf gam rad pos ang

sino = zeros(nb,na);

cangs = cos(angs);
sangs = sin(angs);

% loop over rects
ticker reset
ne = size(rects,1);
for ie = 1:ne
	ticker(mfilename, ie, ne)
	rect = rects(ie,:);

	cx = xscale * rect(1);	wx = rect(3);
	cy = yscale * rect(2);	wy = rect(4);
	eang = deg2rad(rect(5));
	if yscale == -1
		eang = -eang;
	end
	if xscale == -1
		eang = pi - eang;
	end
	val = rect(6);

	% square of projected radius:
	cos_ang_rot = cangs * cos(eang) + sangs * sin(eang);
	sin_ang_rot = sangs * cos(eang) - cangs * sin(eang);
	rp = sqrt((wx * cos_ang_rot).^2 + (wy * sin_ang_rot).^2);

	% radial shift
	sp = cx * cangs + cy * sangs;	
	% sacled distance from center
	dis = (rads - sp) ./ rp;

	% projection angle after affine scaling and rotate
	abs_cos_ang_pi = ones(nb, na);
	nonzeros = sin_ang_rot~=0;
	abs_cos_ang_pi = wx * abs(cos_ang_rot) ./ rp .* nonzeros;

	abs_sin_ang_pi = zeros(nb, na);
	abs_sin_ang_pi = wy * abs(sin_ang_rot) ./ rp .* nonzeros;

	% break points of the trapezoid
	len = 1 ./ max(abs_cos_ang_pi, abs_sin_ang_pi);
	dmax = (abs_cos_ang_pi + abs_sin_ang_pi) / 2;
	dbreak = abs(abs_cos_ang_pi - abs_sin_ang_pi) / 2;
	dmax_break = dmax - dbreak;
	scale = val * wx * wy ./ rp .* len;

	% line integrals at the left part
	leftpart = -dmax < dis & dis < -dbreak;
	sino = sino + scale .* (dis + dmax) ./ dmax_break .* leftpart;
	% line integrals at the middle part
	midpart = abs(dis) <= dbreak;
	sino = sino + scale .* midpart;
	% line integrals at the right part
	rightpart = dbreak < dis & dis < dmax;
	sino = sino + scale .* (dmax - dis) ./ dmax_break .* rightpart;
end


%
% rect_sino_test()
% internal test routine: standard sampling
%
function rect_sino_test

[sys f] = ct_sys;
down = 1;
% only when down=1 and rectangles on the grid,
% wtfmex outperforms DD obviously
ig = image_geom('nx', f.nx, 'ny', f.ny, 'dx', f.pixel_size, ...
	'offset_x', f.center_x, 'offset_y', f.center_y, 'down', down);
rect = [0 0 20*ig.dx 20*ig.dx 0 1];
[xtrue rect]= rect_im(ig, rect, 'oversample', 3);
im pl 2 3
im(1, ig.x, ig.y, xtrue,'xtrue'), cbar

sgf = sino_geom('fan', 'nb', f.nb, 'na', f.na, 'ds', f.ray_spacing, ...
	'dsd', f.dis_src_det, 'dod', f.dis_iso_det, ...
	'offset_s', f.channel_offset, ...
	'strip_width', f.strip_width, ...
	'source_offset', f.source_offset,...
	'dfs', 0, ...%inf, ... % flat fan
	'orbit', f.orbit, 'orbit_start', f.orbit_start, 'down', down);

% analytical sinogram
ya = rect_sino(sgf, rect, 'oversample', 8);
im(2, sgf.s, sgf.ad, ya, 'ya'), cbar

% aspire system 14
% todo: only if user has aspire
% todo use Gtomo2_wtmex
args = arg_pair('system', 14, 'nx', ig.nx, 'ny', ig.ny, ...
	'nb', sgf.nb, 'na', sgf.na, ...
	'support', 'all', ...
	...%	'support', ['file ' f.mask], ...
	'pixel_size', ig.dx, ...
	'orbit', sgf.orbit, 'orbit_start', sgf.orbit_start, ...
	'src_det_dis', sgf.dsd, 'obj2det_x', sgf.dod, ...
	'ray_spacing', sgf.ds, ...
	'channel_offset', sgf.offset_s, ... % fix: -?
	'flip_y', 1, 'scale', 0);

Gdsc = Gtomo2_dscmex(args);

% DD
% todo: only if user has dd
Gdd = Gtomo_dd(sgf, ig, 'nthread', 1);

yd = Gdd * xtrue;
yw = Gdsc * xtrue;

im(3, sgf.s, sgf.ad, yd, 'yd'), cbar,
im(4, sgf.s, sgf.ad, yw, 'yw'), cbar,

t = sprintf('dd err %g%%', max_percent_diff(ya, yd));
im(5, sgf.s, sgf.ad, abs(yd-ya), t), cbar,
t = sprintf('wtfmex err %g%%', max_percent_diff(ya, yw));
im(6, sgf.s, sgf.ad, abs(yw-ya), t), cbar,
nrms(yd, ya)
nrms(yw, ya)
