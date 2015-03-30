function [y,Phi,Phit,Psi,Psit,opt,FBP]=loadCastSim(opt)
    %theta=1:180;     %for phantom
    %theta=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta=1:160;  % Dogandzic2011Asilomar

    if(~isfield(opt,'beamharden')) opt.beamharden=false; end

    opt.trueImg=double(imread('binaryCasting.bmp'));
    conf=ConfigCT();

    daub = 2; dwt_L=6;        %levels of wavelet transform
    opt.rCoeff=[3000 5000 7000 8000 10000 15e3 20000 35000 50000 1e5 5e5]; 
    maskType='CircleMask';
    spark=false;

    conf.PhiMode = 'gpuPrj'; %'parPrj'; %'basic'; %'gpuPrj'; %
    conf.dist = 5000;
    conf.imgSize = size(opt.trueImg,1);
    conf.prjWidth = conf.imgSize;
    conf.prjFull = opt.prjFull;
    conf.prjNum = opt.prjNum;
    conf.dSize = 1;
    conf.effectiveRate = 1;
    conf.Ts = 1;

    [ops.Phi,ops.Phit,ops.FBP]=conf.genOperators();  % without using mask
    if(opt.beamharden)
        symbol={'Fe'};
        densityMap{1}=opt.trueImg;

        [y,args] = genBeamHarden(symbol,densityMap,ops,...
            'showImg',false,'spark', spark);
        opt.trueIota = args.iota(:);
        opt.epsilon = args.epsilon(:);
        opt.trueKappa = args.kappa(:);

        conf.Ts = args.Ts;

        y=-log(y(:)/max(y(:)));
    else
        y = ops.Phi(opt.trueImg);
    end

    y = y(:);
    v = randn(size(y));
    v = v*(norm(y)/sqrt(opt.snr*length(y)));
    y = y + v;

    % remedy for the normalization
    if(opt.beamharden) y=y-min(y); end

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

