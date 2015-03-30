function [y,Phi,Phit,Psi,Psit,opt,FBP]=loadYinyang(opt)
%   n=512;
%   trueImg=zeros(n,n);
%   inFig=rgb2gray(double(imread('yinyang.png')));
%   [~,~,alpha]=imread('yinyang.png');
%   inFig(alpha==0)=0;
%   t1=floor((n-size(inFig,1))/2);
%   t2=ceil((n-size(inFig,1))/2);
%   trueImg(t1+1:n-t2,t1+1:n-t2)=inFig;

    if(~isfield(opt,'beamharden')) opt.beamharden=false; end

    trueImg=load('yinyang1.mat'); opt.trueImg=trueImg.trueImg;
    conf=ConfigCT();

    daub = 2; dwt_L=6;        %levels of wavelet transform
    maskType='CircleMask';
    spark=false;

    conf.PhiMode = 'gpuPrj'; %'parPrj'; %'basic'; %'gpuPrj'; %
    conf.dist = 2000;
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

        [y,args] = genBeamHarden('showImg',false,...
            'spark', spark, 'densityMap',densityMap, ...
            'operators',ops, 'symbol',symbol);
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
        mask = Utils.getCircularMask(n);
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',n,dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,dwt_L,daub,wvltName);
        end
    end
    opt.mask = mask;
    opt.trueAlpha=opt.trueImg(maskIdx);

    [Phi,Phit,FBP]=conf.genOperators(mask);
    [Psi,Psit]=Utils.getPsiPsit(daub,dwt_L,mask,maskk);

    fprintf('Configuration Finished!\n');
end

