function loadCastSim(obj)
    %theta=1:180;     %for phantom
    %theta=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta=1:160;  % Dogandzic2011Asilomar

    obj.trueImg=double(imread('binaryCasting.bmp'));
    [obj.CTdata,args] = genBeamHarden('showImg',false,...
        'spark', obj.spark, 'trueImg',obj.trueImg, ...
        'prjNum', obj.prjNum, 'prjFull', obj.prjFull,...
        'PhiMode',obj.PhiModeGen);
    obj.trueIota = args.iota(:);
    obj.epsilon = args.epsilon(:);
    obj.trueKappa = args.kappa(:);
    obj.Ts = args.Ts;

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
end

