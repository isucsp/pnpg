% This configuration file is used for RealCT data and the analytically 
% simulated Phantom data. The sinogram is first transformed to frequency 
% domain, if the length of every projection is N, then the reconstructed image 
% should have a size of NxN.

% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.2 $ $Date: Wed 05 Feb 2014 05:09:42 PM CST
% v_0.2:        change the structure to class for easy control;

classdef ConfigCT < handle
    properties
        imageName
        maskType = 'CircleMask';
        operType = 'nufft'; %'gpu'; %'cpu'

        % for generating polychromatic measurements
        spark   

        % generated data for further usage
        trueImg
        Ts

    end 
    methods
        function genOperators(obj)
            switch lower(obj.operType)
                case 'nufft'
                case 'gpu'
                case 'cpu'
            end
        end
        function gen
Num_proj=size(CTdata,2);
theta_full=linspace(0,pi-pi/Num_proj,Num_proj); %The angle of projections
theta=theta_full(theta_idx);
CTdata=CTdata(:,theta_idx);
[Num_pixel,Num_proj]=size(CTdata);

temp=2^(ceil(log2(Num_pixel)));
CTdata=[zeros(ceil((temp-Num_pixel)/2),Num_proj); CTdata;...
    zeros(floor((temp-Num_pixel)/2),Num_proj)];
[Num_pixel,Num_proj]=size(CTdata);

mx=Num_pixel;
my=mx; m=mx*my;
m_2D=[mx, my];
J=[1,1]*3;                       % NUFFT interpolation neighborhood
K=2.^ceil(log2(m_2D*2));         % oversampling rate

r=pi*linspace(-1,1-2/Num_pixel,Num_pixel)';

xc=r*cos(theta(:)');
yc=r*sin(theta(:)');
om=[yc(:), xc(:)];
st=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');
st.Num_pixel=Num_pixel;
st.Num_proj=Num_proj;

y=CTdata(:);
N=length(y);
Imea=exp(-y);

% Zero freq at f_coeff(Num_pixel/2+1)
f_coeff=designFilter('renliang1',Num_pixel,Ts);

Mask=double(Mask~=0);
%figure; showImg(Mask);
wvltIdx=find(maskk~=0);
p_I=length(wvltIdx);
maskIdx=find(Mask~=0);
p_M=length(maskIdx);

fprintf('Generating Func handles...\n');

%Sampling operator
W=@(z) midwt(z,wav,dwt_L);
Wt=@(z) mdwt(z,wav,dwt_L);

Psi=@(s) PsiFunc(s,W);
PsiM=@(s) PsiFunc(s,W,m,maskIdx,wvltIdx);
Psit=@(s) PsitFunc(s,Wt);
PsiMt=@(s) PsitFunc(s,Wt,m,maskIdx,wvltIdx);

switch lower(PhiPhitMode)
    case 'basic'
        % The most basic methods without anything here
        switch 0
            case 0
                Phi=@(s) PhiFunc51(s,f_coeff,st,mx,Ts);
                PhiM=@(s) PhiFunc51(s,f_coeff,st,mx,Ts,maskIdx);
                Phit=@(s) PhitFunc51(s,f_coeff,st,mx,Ts);
                PhiMt=@(s) PhitFunc51(s,f_coeff,st,mx,Ts,maskIdx);
                FBP=@(s) FBPFunc6(s,theta_full,theta_idx,Ts);
            case 1
                conf.bw=1; conf.nc=mx; conf.nr=my; conf.prjWidth=Num_pixel;
                conf.theta=theta*180/pi;
                temp=1:numel(Img2D);
                Phi= @(s) mParPrj(s,temp-1,conf,'forward')*Ts;
                Phit =@(s) mParPrj(s,temp-1,conf,'backward')*Ts;
                PhiM=@(s) mParPrj(s,maskIdx-1,conf,'forward')*Ts;
                PhiMt=@(s) mParPrj(s,maskIdx-1,conf,'backward')*Ts;
                FBP=@(s) FBPFunc7(s,theta_full,theta_idx,Ts,maskIdx)*Ts;
            case 2
                conf.n=mx; conf.prjWidth=Num_pixel; conf.np=Num_proj; 
                conf.prjFull=360; conf.dSize=(mx-1)/(Num_pixel+1); conf.effectiveRate=1; 
                conf.d=dist;

                mPrj(0,conf,'config');
                mPrj(0,0,'showConf');

                Phi= @(s) mPrj(reshape(single(s),conf.n,[])',conf,'forward')*Ts;
                Phit =@(s) reshape(reshape(mPrj(single(s(:)),conf,'backward'),conf.n,[])'*Ts,[],1);
                PhiM=@(s) mPrj(single(maskFunc(s,maskIdx,conf.n))',conf,'forward')*Ts;
                PhiMt=@(s) maskFunc(reshape(mPrj(single(s(:)),conf,'backward'),conf.n,[])'*Ts,maskIdx);
                FBP=@(s) FBPFunc8(s,conf,Ts,maskIdx)*Ts;
        end

    case 'weighted'
        % Fessler's weighted methods
        weight=exp(-y);
        %weight=exp(-ones(size(y)));
        weight=sqrt(weight/max(weight(:)));
        y=y.*weight;

        Phi=@(s) PhiFunc51(s,f_coeff,st,mx,Ts).*weight(:);
        PhiM=@(s) PhiFunc51(s,f_coeff,st,mx,Ts,maskIdx).*weight(:);
        Phit=@(s) PhitFunc51(s.*weight(:),f_coeff,st,mx,Ts);
        PhiMt=@(s) PhitFunc51(s.*weight(:),f_coeff,st,mx,...
            Ts,maskIdx);

        FBP=@(s) FBPFunc6(s./weight,theta_idx,Ts);

    case 'filtered'
        % Sqrt filtered methods
        y=reshape(y(:),Num_pixel,Num_proj);
        y=[zeros(Num_pixel/2,Num_proj); y; zeros(Num_pixel/2,Num_proj)];
        y=fft(fftshift(y,1))*Ts;
        y=y.*repmat(sqrt(f_coeff),1,Num_proj);
        y=fftshift(ifft(y),1)/Ts;
        y=y(Num_pixel/2+1:Num_pixel/2+Num_pixel,:);
        y=real(y(:));

        Phi=@(s) PhiFunc2(s,f_coeff,stFwd,Num_pixel,Ts);
        PhiM=@(s) PhiFunc2(s,f_coeff,stFwd,Num_pixel,Ts,maskIdx);
        Phit=@(s) PhitFunc2(s,f_coeff,stFwd,Num_pixel,Ts);
        PhiMt=@(s) PhitFunc2(s,f_coeff,stFwd,Num_pixel,Ts,maskIdx);

        FBP=Phit;

    otherwise
        fprintf('Wrong mode for PhiPhitMode: %s\n',PhiPhitMode);
        return;
end

H=@(s) Phi(Psi(s));
HM=@(s) PhiM(PsiM(s));
Ht=@(s) Psit(Phit(s));
HMt=@(s) PsiMt(PhiMt(s));

scale=max(Img2D(:))-min(Img2D(:));
scaleM=max(Img2D(maskIdx))-min(Img2D(maskIdx));
normy=norm(y);

fprintf('Configuration Finished!\n');

if(0)
    testTranspose(Phi,Phit,N,m,'Phi');
    testTranspose(PhiM,PhiMt,N,p_M,'PhiM');
%   testTranspose(Psi,Psit,m,m,'Psi');
%   testTranspose(PsiM,PsiMt,p_M,p_I,'PsiM');
end

c=8.682362e-03;
mc=1;
mc=0.7;
mc=6.8195e-03;    % For full projection
%mc=1.307885e+01;
if(0)
    %c=expectHHt(H,Ht,N,m,'H');
    c=expectHHt(Phi,Phit,N,m,'Phi');
    %mc=expectHHt(HM,HMt,N,p_I,'HM');
    mc=expectHHt(PhiM,PhiMt,N,p_M,'PhiM');
end

img=zeros(size(Img2D));
img(567:570,787:790)=1;
%y=Phi(img);

        function loadMeasurements(obj)
            fprintf('Loading data...\n');
            switch lower(obj.imageName)
                case 'phantom'
                    loadPhantom(obj);
                case lower('castSim')
                    loadCastSim(obj);
                case lower('twoMaterials')
                    loadTwoMaterials(obj);
                case 'realct'
                    loadRealCT(obj);
                case 'pellet'
                    loadPellet(obj);
            end
        end
        function loadCastSim(obj)
            obj.Ts=0.008;
            theta_idx=1:180;     %for phantom
            %theta_idx=[1:10, 21:100, 111:180]; % Kun2012TSP cut
            %theta_idx=1:160;  % Dogandzic2011Asilomar

            obj.trueImg=double(imread('binaryCasting.bmp'));
            obj.CTdata = genBeamHarden('showImg',false, 'spark', obj.spark,...
                'trueImg',obj.trueImg);
            obj.CTdata=-log(obj.CTdata/max(obj.CTdata(:)));

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
            obj.wav=daubcqf(2);
            obj.dwt_L=6;        %levels of wavelet transform
            obj.rCoeff=[3000 5000 7000 8000 10000 15e3 20000 35000 50000 1e5 5e5]; 
        end
        function loadTwoMaterials()
            Ts=0.008;
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

        function fan2parallel()
            [temp1, ploc1, ptheta1]=...
                fan2paraM(CTdata,dist*Ts,'FanSensorGeometry','line',...
                'FanSensorSpacing',Ts,'ParallelCoverage','halfcycle',...
                'Interpolation','pchip','ParallelRotationIncrement',paraInc,...
                'PerpenCenter',perpenCenter,'RotationCenter',rotationCenter,...
                'ParallelSensorSpacing',Ts*ds);

            [temp2, ploc2, ptheta2]=...
                fan2paraM(CTdata(:,[181:360 1:180]),dist*Ts,...
                'FanSensorGeometry','line',...
                'FanSensorSpacing',Ts,'ParallelCoverage','halfcycle',...
                'Interpolation','pchip','ParallelRotationIncrement',paraInc,...
                'PerpenCenter',perpenCenter,'RotationCenter',rotationCenter,...
                'ParallelSensorSpacing',Ts*ds); %fss*ds);  %

            temp2=temp2(end:-1:1,:);
            ploc2=-ploc2(end:-1:1);
            lowerBound=max(ploc1(1),ploc2(1));
            upperBound=min(ploc1(end),ploc2(end));
            if(lowerBound==ploc1(1))
                ploc=ploc2; ploc2=ploc1;
                CTdata=temp2; temp2=temp1;
            else
                ploc=ploc1;
                CTdata=temp1;
            end
            idx1=find(ploc<lowerBound);
            idx2=find((ploc>=lowerBound) & (ploc <=upperBound));
            idx3=idx1+length(ploc);
            ploc(idx3)=ploc2(end-length(idx3)+1:end);
            CTdata=[CTdata; zeros(length(idx3),size(CTdata,2))];
            CTdata([idx2;idx3],:)=CTdata([idx2;idx3],:)+temp2;
            CTdata(idx2,:)=CTdata(idx2,:)/2;
            % Try to reduce the affect of beam harden
            % filt = (.74 + .26 * cos(pi*(-Num_pixel/2:Num_pixel/2-1)/Num_pixel));
            % CTdata=CTdata.*repmat(filt(:),1,Num_proj);
            % aaa=10^-0.4;
            % CTdata=(exp(aaa*CTdata)-1)/aaa;
        end
    end
end

