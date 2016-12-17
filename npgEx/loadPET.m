function [y,Phif,Phitf,Psi,Psit,fbpfunc,opt]=loadPET(totalCnt,opt,seed)
    if(~exist('seed','var')) seed=0; end
    if(~exist('totalCnt','var')) totalCnt=1e7; end
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',seed));

    % true emission image
    ig = image_geom('nx', 128, 'ny', 128, 'fov', 500);
    % emission image
    xtrue = read_zubal_emis('nx', ig.nx, 'ny', ig.ny);
    % attenuation map
    mumap = read_zubal_attn('nx', ig.nx, 'ny', ig.ny);

    x=xtrue;
    if(~exist('opt','var'))
        [~,l1norm]=Utils.testSparsity(x);
        opt.l1norm=l1norm;
        opt.x=x;
        figure; showImg(x);
        return;
    end

    dwt_L=6; daub=6;

    if(isfield(opt,'mask') && isempty(opt.mask))
        mask=[]; maskk=[];
    else
        % reconstruction mask (which pixels do we estimate?)
        mask = Utils.getCircularMask(ig.nx);
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',ig.nx,dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,dwt_L,daub,wvltName);
        end
        opt.mask=mask; opt.maskk=maskk;
        ig.mask = mask > 0;
        % ig.mask = ig.circ(240, 240) > 0;
        % showImg(ig.mask + xtrue); title('support mask + xtrue');
    end

    % system matrix G
    sg = sino_geom('par', 'nb', ig.nx, 'na', ig.ny*0+90, ...
        'dr', 500/(ig.nx));

    % simple strip-integral system model
    G = Gtomo2_strip(sg, ig);

    % noisy measurements
    proj = G * xtrue;
    li = G * mumap;
    % printm('Maximum line integral = %g', max(li(:)))
    if ~isvar('f.count'), f.count = totalCnt; end
    % detector efficiency variations per CTI 931 PET scanner
    ci = exp(0.3 * randn(size(proj)));
    ci = ci .* exp(-li);
    ci = f.count / sum(ci(:) .* proj(:)) * ci;
    ytrue = ci .* proj;
    if ~isvar('f.randpercent')
        f.randpercent = 10;
    end
    ri = f.randpercent / 100 * mean(ytrue(:)) * sg.ones;
    if(isfield(opt,'noBackground') && opt.noBackground)
        ri=ri*0;
    end
    yi = poisson(ytrue + ri);

    y=yi;

    % FBP reconstruction
    xfbp = em_fbp(sg, ig, yi, ci, ri);
    fbpfunc = @(yyy) double(em_fbp(sg, ig, reshape(yyy,size(ci)), ci, ri));

    Phif = @(xx) ci.*(G*xx);
    Phitf = @(xx) G'*(ci.*xx);

    % caution: don't use the wavelet tools from matlab, it is slow
    wav=daubcqf(daub);
    W=@(z) midwt(z,wav,dwt_L);
    Wt=@(z) mdwt(z,wav,dwt_L);

    if(~isempty(mask))
        maskIdx = find(mask~=0);
        wvltIdx = find(maskk~=0);
        Psi = @(s) maskFunc(W (maskFunc(s,wvltIdx,ig.nx)),maskIdx);
        Psit= @(x) maskFunc(Wt(maskFunc(x,maskIdx,ig.nx)),wvltIdx);
        opt.trueX=xtrue(maskIdx);
    else
        Psi = W;
        Psit= Wt;
        opt.trueX=xtrue;
    end

    opt.bb=ri;
    % opt.ci=ci(:);
    % opt.G=G;
end

