function loadGlassBeadsSim(obj,snr,noiseType)
    %theta=1:180;     %for phantom
    %theta=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta=1:160;  % Dogandzic2011Asilomar

    if(nargin<=2) noiseType='poissonLogLink'; end
    obj.trueImg=double(imread('glassbeads.bmp'));
    obj.trueImg = obj.trueImg/max(obj.trueImg(:));

    daub = 2;
    obj.wav=daubcqf(daub);
    obj.dwt_L=4;        %levels of wavelet transform
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

                obj.Ts=1;
                genOperators(obj, obj.PhiModeGen);
                % suppose Φx \in [a,b], we want to map I_0 exp(-Φx) to [A,B]
                y = obj.Phi(obj.trueImg(mask~=0));
                a=min(y); b=max(y); A=10; B=2^12;
                scale=(log(B)-log(A))/(b-a);
                obj.I0=exp( (b*log(B) - a*log(A))/(b-a) );

                obj.Ts=obj.Ts*scale;
                y = obj.I0*exp(-y*scale);
                obj.y = poissrnd(y);
            case 'gaussian'
                obj.Ts=1;
                genOperators(obj, obj.PhiModeGen);
                obj.y = obj.Phi(obj.trueImg(mask~=0));
        end
    end
end

