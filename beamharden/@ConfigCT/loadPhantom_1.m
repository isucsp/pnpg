function loadPhantom_1(obj)
    %theta=1:180;     %for phantom
    %theta=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta=1:160;  % Dogandzic2011Asilomar

    obj.trueImg=zeros(obj.imgSize,obj.imgSize);
    [x,y] = meshgrid(1:obj.imgSize);
    idx = sqrt((x-0.52*obj.imgSize).^2+(y-obj.imgSize*0.52).^2) < 0.32*obj.imgSize;
    obj.trueImg(idx)=1;
    idx = sqrt((x-0.62*obj.imgSize).^2+(y-obj.imgSize*0.62).^2) < 0.12*obj.imgSize;
    obj.trueImg(idx)=0;

    daub = 2;
    obj.wav=daubcqf(daub);
    obj.dwt_L=6;        %levels of wavelet transform
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

    if(obj.beamharden)
        [obj.CTdata,args] = genBeamHarden('showImg',false,...
            'spark', obj.spark, 'trueImg',obj.trueImg, ...
            'prjNum', obj.prjNum, 'prjFull', obj.prjFull,...
            'PhiMode',obj.PhiModeGen);
        obj.trueIota = args.iota(:);
        obj.epsilon = args.epsilon(:);
        obj.trueKappa = args.kappa(:);
        obj.Ts = args.Ts;
        
        % before this, make sure max(CTdata)==1
        temp=size(obj.CTdata,1);
        obj.CTdata=[zeros(ceil((obj.prjWidth-temp)/2),obj.prjNum);...
            obj.CTdata;...
            zeros(floor((obj.prjWidth-temp)/2),obj.prjNum)];
        obj.y=-log(obj.CTdata(:)/max(obj.CTdata(:)));
    else
        obj.Ts=1;
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

