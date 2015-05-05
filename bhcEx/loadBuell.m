function [y,Phi,Phit,Psi,Psit,opt,FBP]=loadBuell(opt)

    y=load('buell_1a225.mat');
    y=y.data;

    M=size(y,2);
    y=y(:,1:M/opt.prjFull:M/opt.prjFull*opt.prjNum);

    conf=ConfigCT();

    daub = 2; dwt_L=9;        %levels of wavelet transform
    maskType='none';

    conf.PhiMode = 'gpuPrj'; %'parPrj'; %'basic'; %'gpuPrj'; %
    conf.dist = 7765;
    conf.imgSize = size(y,1);
    conf.prjWidth = conf.imgSize;
    conf.prjFull = opt.prjFull;
    conf.prjNum = opt.prjNum;
    conf.dSize = 1;
    conf.effectiveRate = 1;
    conf.Ts =1e-2;

    detectorBitWidth=16;

    if(strcmp(maskType,'CircleMask'))
        % reconstruction mask (which pixels do we estimate?)
        mask = Utils.getCircularMask(conf.imgSize);
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',conf.imgSize,dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,dwt_L,daub,wvltName);
        end
    else
        mask = []; maskk=[];
    end
    [Phi,Phit,FBP]=conf.genOperators(mask);
    [Psi,Psit]=Utils.getPsiPsit(daub,dwt_L,mask,maskk);

    fprintf('Configuration Finished!\n');
end

