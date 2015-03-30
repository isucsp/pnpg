function loadPellet(obj)
    Ts=0.01*0.8;

    theta_idx=1:180;
    rCoeff=[2000 3000 5000 7000 8000 10000 20000 35000 50000 1e5 5e5]; %1024x1024RealCT
    %Pre-processing data
    % For *non-calibrated* data
    %rotationCenter=1008.7; dist=7248.2; perpenCenter=695;
    %load 'fuelPellet/test036_0223.mat'
    %ds=1.96;

    % For *calibrated* data
    rotationCenter=883.6; dist=8855; perpenCenter=610;
    % FAN_BEAM_PARAMETERS       1086.6030       24.4266
    % COR_PARAMETERS            -14735.476562   0.001953
    load '../fuelPellet.db/test036_Calibrated_0302.mat'
    ds=1.7214;
    %load(['/home/renliang/research/fuelPellet.db/' ...
    %    'test036_Calibrated_0150.mat']);
    wrap=1; %50; %

    wav=daubcqf(6);
    dwt_L=8;        %levels of wavelet transform
    %       wav=daubcqf(2);
    %       dwt_L=6;        %levels of wavelet transform

    load 'RealCT.mat';

    fprintf('Loading %s...\n',obj.maskType);
    if(strcmp(obj.maskType,'/CircleMask'))
        load('Mask1024Circle.mat');
        load('MaskWvlt1024CircleL8D6.mat');
    else 
        load(sprintf('RealCTMask_%02d.mat',maskNo));
        load(sprintf('RealCTMaskWvlt_%02d.mat',maskNo));
    end

    CTdata_mtrx=CTdata;
    CTdata=-log(CTdata/max(CTdata(:)));
    CTdata=CTdata(:,end:-1:1);
    CTdata=CTdata(:,[wrap:360 1:wrap-1]);

    paraInc=1;

    if change
        fan2parallel();
    end
end

