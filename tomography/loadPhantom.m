function loadPhantom(obj,snr,noiseType)
    % use 'Modified Shepp-Logain' Phantom for the truth
    obj.trueImg=phantom(obj.imgSize);
    obj.trueImg(obj.trueImg<0)=0;
    obj.prjWidth= obj.imgSize/obj.dSize;

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
        error('Phantom with beamhardening is not implementaed');
    else
        switch lower(noiseType)
            case lower('poissonLogLink')
                % u = 2e-4 is good choice
                scale = 1e-1;           % determins Ts
                obj.Ts=1;
                genOperators(obj, obj.PhiModeGen);

                PhiAlpha = obj.Phi(obj.trueImg(mask~=0));
                obj.Ts = log(1/scale)/max(PhiAlpha);

                PhiAlpha = PhiAlpha*obj.Ts;
                y = exp(-PhiAlpha);
                I0 = (snr-1)*sum(y)/sum(y.^2);
                obj.y = poissrnd(I0*y);
                obj.y(obj.y==0) = 1;
            case 'poisson'
                obj.Ts=1;
                genOperators(obj, obj.PhiModeGen);

                PhiAlpha=obj.Phi(obj.trueImg(mask~=0));
                a = snr*sum(PhiAlpha)/sum(PhiAlpha.^2);
                obj.trueImg = obj.trueImg*a;
                obj.y = poissrnd(a*PhiAlpha);
            case 'gaussian'
                obj.Ts=1;
                genOperators(obj, obj.PhiModeGen);
                obj.y = obj.Phi(obj.trueImg(mask~=0));

                v = randn(size(obj.y));
                v = v*(norm(obj.y)/sqrt(snr*length(obj.y)));
                obj.y = obj.y + v;
        end
    end
end

