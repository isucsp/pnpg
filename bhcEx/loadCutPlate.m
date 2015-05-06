function [y,Phi,Phit,Psi,Psit,opt,FBP]=loadBuell(opt)

    y=load('cut_plate_a_220.mat');
    y=y.data;
    y=y-min(y(:));
    if(mod(size(y,1),2)==1) y=[zeros(1,size(y,2)); y];

    M=size(y,2);
    y=y(:,1:M/opt.prjFull:M/opt.prjFull*opt.prjNum);

    conf=ConfigCT();

    daub = 2; dwt_L=6;        %levels of wavelet transform
    maskType='CircleMask';

    conf.PhiMode = 'gpuPrj'; %'parPrj'; %'basic'; %'gpuPrj'; %
    conf.imgSize = 1024;
    conf.prjWidth = size(y,1);
    conf.prjFull = opt.prjFull;
    conf.prjNum = opt.prjNum;
    conf.dSize = conf.imgSize/conf.prjWidth;
    conf.effectiveRate = 1;
    conf.dist = 8900/conf.prjWidth*conf.imgSize;
    conf.Ts =1e-2;

    if(strcmpi(maskType,'CircleMask'))
        % reconstruction mask (which pixels do we estimate?)
        mask = Utils.getCircularMask(conf.imgSize);
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',conf.imgSize,dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,dwt_L,daub,wvltName);
        end
    elseif(strcmpi(maskType,'none'))
        mask = []; maskk=[];
    end
    [Phi,Phit,FBP]=conf.genOperators(mask);
    [Psi,Psit]=Utils.getPsiPsit(daub,dwt_L,mask,maskk);
    opt.mask=mask; opt.maskk=maskk;

    fprintf('Configuration Finished!\n');
end

