function loadCastSim(obj)
    obj.Ts=0.008;
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

    if(strcmp(obj.maskType,'CircleMask'))
        load('Mask1024Circle.mat');
        load('MaskWvlt1024CircleL8D6.mat');
    elseif(strcmp(obj.maskType,'cvxHull'))
        load(sprintf('RealCTMask_%02d.mat',2));
        load(sprintf('RealCTMaskWvlt_%02d.mat',2));
    else
        load(sprintf('RealCTMask_%02d.mat',3));
        load(sprintf('RealCTMaskWvlt_%02d.mat',3));
    end
    obj.mask = Mask; obj.maskk= maskk;
    obj.wav=daubcqf(2);
    obj.dwt_L=6;        %levels of wavelet transform
    obj.rCoeff=[3000 5000 7000 8000 10000 15e3 20000 35000 50000 1e5 5e5]; 
end

