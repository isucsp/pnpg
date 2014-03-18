function loadPhantom(obj)
    % use 'Modified Shepp-Logain' Phantom for the truth
    obj.trueImg=phantom(obj.imgSize);

    %levels of wavelet transform
    if(obj.imgSize<=512)
        obj.dwt_L=6;    %for  512x512  Phantom image
        obj.rCoeff=[6e3 7e3 8e3 1e4 5e4 2e5];   %512x512Phantom
    else
        obj.rCoeff=[15000 16000 17000 18e3 19e3 2e4 3e4 5e4 4e4 6e4];
        %1024x1024Phantom
        obj.dwt_L=6;   %for 1024x1024 Phantom image
    end
    daub = 2;
    obj.wav=daubcqf(daub);

    if(strcmp(obj.maskType,'CircleMask'))
        mask = zeros(obj.imgSize);
        [x,y]=meshgrid(0:obj.imgSize-1);
        idx = sqrt((x-obj.imgSize/2).^2+(y-obj.imgSize/2).^2)<(obj.imgSize/2-1);
        mask(idx)=1;
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',obj.imgSize,obj.dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,obj.dwt_L,daub,wvltName);
        end
    end
    obj.mask = mask; obj.maskk= maskk;

    if(obj.beamharden)
        error('\nPhantom with beamhardening is not implementaed\n');
    else
        obj.Ts=0.1;
        conf.n=obj.imgSize; conf.prjWidth=obj.prjWidth;
        conf.np=obj.prjNum; conf.prjFull=obj.prjFull;
        conf.dSize=obj.dSize; %(n-1)/(Num_pixel+1);
        conf.effectiveRate=obj.effectiveRate;
        conf.d=obj.dist;

        mPrj(0,conf,'config');
        maskIdx = find(obj.mask~=0);
        Phi =@(s) mPrj(maskFunc(s,maskIdx,conf.n),0,'forward')*obj.Ts;
        obj.y = Phi(obj.trueImg(obj.mask~=0));
    end
end
