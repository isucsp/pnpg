function loadRealCT(obj)
    Ts=0.01*0.8;

    theta_idx=1:180;  %[1,79]; %
    %theta_idx=[1:10, 21:100, 111:180]; % Kun2012TSP cut
    %theta_idx=1:160;  % Dogandzic2011Asilomar

    %Pre-processing data
    center_loc=691; %692;
    dist=8597; %8797; %
    %FAN_BEAM_PARAMETERS 1098.8800	    536.0000
    %COR_PARAMETERS -28.614487	    0.001389
    load 'CTdata_220.mat'
    %1024x1024RealCT
    rCoeff=[3000 5000 7000 8000 10000 15e3 20000 35000 50000 1e5 5e5]; 
    % In asilomar paper: DORE 2e4 maskDORE 15e3
    % In Kun's Juarnal: DORE 5000 and 2e4
    wrap=1;

    wav=daubcqf(6);
    dwt_L=8;        %levels of wavelet transform

    load 'RealCT.mat';

    fprintf('Loading %s...\n',obj.maskType);
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
    CTdata(CTdata==max(CTdata(:)))=max(CTdata(:))*1.0;
    CTdata=-log(CTdata/max(CTdata(:)));
    CTdata=CTdata(:,end:-1:1);
    CTdata=CTdata(:,[wrap:360 1:wrap-1]);
    CTdata=CTdata(81:2*center_loc-81,:);
    %CTdata(end-98:end,:)=[];
    %CTdata(1:98,:)=[];

    ds=1.1924;
    paraInc=1;
    %           CTdata=CTdata_mtrx;
    perpenCenter=(size(CTdata,1)+1)/2;
    rotationCenter=perpenCenter;
    %           ds=0.8419;
end
