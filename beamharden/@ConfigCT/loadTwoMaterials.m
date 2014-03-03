function loadTwoMaterials(obj)
    obj.Ts=0.008;
    theta_idx=1:180;     %for phantom
    %theta_idx=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta_idx=1:160;  % Dogandzic2011Asilomar

    if(spark)
        load 'CTdata_220TwoMaterialsSpark.mat'
    else
        load 'CTdata_220TwoMaterials.mat'
    end
    load('twoMaterialCasting.mat'); %load Img2D
    if(strcmp(obj.maskType,'/CircleMask'))
        load('Mask1024Circle.mat');
        load('MaskWvlt1024CircleL8D6.mat');
    elseif(strcmp(obj.maskType,'/cvxHull'))
        load(sprintf('RealCTMask_%02d.mat',2));
        load(sprintf('RealCTMaskWvlt_%02d.mat',2));
    else
        load(sprintf('RealCTMask_%02d.mat',3));
        load(sprintf('RealCTMaskWvlt_%02d.mat',3));
    end

    CTdata=CTdata_mtrx;
    CTdata=-log(CTdata/max(CTdata(:)));

    wav=daubcqf(6);
    dwt_L=8;        %levels of wavelet transform
    rCoeff=[3000 5000 7000 8000 10000 15e3 20000 35000 50000 1e5 5e5]; 
end
