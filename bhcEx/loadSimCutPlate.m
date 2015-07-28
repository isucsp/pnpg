function [y,Phi,Phit,Psi,Psit,opt,FBP,mask]=loadSimCutPlate(opt)
    RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
    if(~isfield(opt,'beamharden')) opt.beamharden=false; end

    trueImg=double(imread('simCutPlate.png'));
    trueImg=trueImg/max(trueImg(:));
    opt.trueImg=trueImg;

    conf=ConfigCT();

    daub = 2; dwt_L=6;        %levels of wavelet transform
    maskType='CircleMask';

    conf.PhiMode = 'gpuPrj'; %'parPrj'; %'basic'; %'gpuPrj'; %
    conf.imgSize = 1024;
    conf.prjWidth = 1370;
    conf.prjFull = opt.prjFull;
    conf.prjNum = opt.prjNum;
    conf.dSize = conf.imgSize/conf.prjWidth;
    conf.effectiveRate = 1;
    conf.dist = 8900/conf.prjWidth*conf.imgSize;
    conf.Ts =1e-2;

    detectorBitWidth=16;

    [ops.Phi,ops.Phit,ops.FBP]=conf.genOperators();  % without using mask
    if(opt.beamharden)
        symbol={'Fe'};
        densityMap{1}=opt.trueImg;

        [y,args] = genBeamHarden(symbol,densityMap,ops,...
            'showImg',false);
        opt.iota = args.iota(:);
        opt.epsilon = args.epsilon(:);
        opt.kappa = args.kappa(:);

        opt.trueImg=opt.trueImg*args.density;
        conf.Ts = args.Ts;

        y=-log(y(:)/max(y(:)));

        %  Poisson measurements
        Imea = 2^detectorBitWidth * exp(-y);
        Imea = poissrnd(Imea);

        %%%  make the data to be saturated.
        if(isfield(opt,'saturated') && opt.saturated)
            Imea=min(Imea,max(Imea(:))*0.8);
        end

        y=-log(Imea/max(Imea));
    else
        y = ops.Phi(opt.trueImg);

        % use Gaussian noise
        v = randn(size(y));
        v = v*(norm(y)/sqrt(opt.snr*length(y)));
        y = y + v;

        % remedy for the normalization, use only for log link
        % if(opt.beamharden) y=y-min(y); end
    end

    if(strcmp(maskType,'CircleMask'))
        % reconstruction mask (which pixels do we estimate?)
        mask = Utils.getCircularMask(conf.imgSize);
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',conf.imgSize,dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,dwt_L,daub,wvltName);
        end
    end
    opt.mask = mask;
    opt.trueAlpha=opt.trueImg(mask~=0);

    [Phi,Phit,FBP]=conf.genOperators(mask);
    [Psi,Psit]=Utils.getPsiPsit(daub,dwt_L,mask,maskk);

    fprintf('Configuration Finished!\n');
end

