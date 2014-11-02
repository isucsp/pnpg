function loadWrist(obj,snr,noiseType)
    %theta=1:180;     %for phantom
    %theta=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta=1:160;  % Dogandzic2011Asilomar

    load('wrist.mat');
    obj.trueImg=Img2D;
    obj.imgSize=sqrt(length(obj.trueImg(:)));
    keyboard

    daub = 6;
    obj.wav=daubcqf(daub);
    obj.dwt_L=3;        %levels of wavelet transform
    obj.rCoeff=[3000 5000 7000 8000 10000 15e3 20000 35000 50000 1e5 5e5]; 

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
    
    if(~obj.beamharden)
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

                PhiAlpha=obj.Phi(obj.trueImg);
                a = opt.snr*sum(y)/sum(y.^2);
                obj.Ts = obj.Ts*a;

                obj.y = poissrnd(a*PhiAlpha);
            case 'gaussian'
                obj.Ts=1;
                genOperators(obj, obj.PhiModeGen);
                obj.y = obj.Phi(obj.trueImg(mask~=0));
        end
    end
end

