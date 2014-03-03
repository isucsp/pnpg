function loadPhantom(obj)
    Ts=0.1;
    theta_idx=1:4:720;     %for phantom
    theta_idx=theta_idx(1:180);
    load 'PhantomPb511.mat';
    load 'Phantom512.mat';
    if(strcmp(obj.maskType,'/CircleMask'))
        load 'Mask512Circle.mat'
        load 'MaskWvlt512Circle.mat'
    else 
        load 'PhantomMask512.mat';
        load 'PhantomMaskWvlt512.mat' ;
    end
    %levels of wavelet transform
    if(size(Img2D,1)==512)
        dwt_L=5;    %for  512x512  Phantom image
        rCoeff=[6e3 7e3 8e3 1e4 5e4 2e5];   %512x512Phantom
    else
        rCoeff=[15000 16000 17000 18e3 19e3 2e4 3e4 5e4 4e4 6e4];
        %1024x1024Phantom
        dwt_L=6;   %for 1024x1024 Phantom image
    end
    wav=daubcqf(2);
end
