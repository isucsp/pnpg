% This configuration file is used for RealCT data and the analytically 
% simulated Phantom data. The sinogram is first transformed to frequency 
% domain, if the length of every projection is N, then the reconstructed image 
% should have a size of NxN.

% Author: Renliang Gu (renliang@iastate.edu)
% $Revision: 0.2 $ $Date: Fri 21 Feb 2014 09:24:14 PM CST
% v_0.2:        change the structure to class for easy control;

classdef ConfigCT < handle
    properties
        imageName = 'castSim'; %'phantom' %'twoMaterials'; %'realct'; %'pellet'; %
        maskType = 'CircleMask';

        PhiMode = 'basic'; %'gpuPrj'; %'cpuPrj'; %'parPrj';
        PhiModeGen = 'parPrj';
        imgSize = 1024;
        prjWidth = 1024;
        prjNum = 180;
        prjFull = 360;
        dSize = 1;
        effectiveRate = 1;
        dist = 0;       % default to be parallel beam

        Phi
        Phit
        Psi
        Psit
        FBP

        % for wavelet configuration
        wav = daubcqf(2);
        dwt_L=6;        %levels of wavelet transform

        rCoeff

        CTdata      % the raw data from the measurement in the matrix form
        y
        mask
        maskk

        % for generating polychromatic measurements
        trueIota
        epsilon
        trueKappa
        spark   

        % generated data for further usage
        trueImg
        Ts

        % parameters for operators

    end 
    methods
        function obj = ConfigCT(in,mt,ot)
            if(nargin>0)
                obj.imageName = in;
                if(nargin>1)
                    obj.maskType = mt;
                    if(nargin>2)
                        obj.opTye = ot;
                    end
                end
            end
        end

        function opt=setup(obj,opt)
            loadMeasurements(obj);
            genOperators(obj);

            opt.trueKappa= obj.trueKappa;
            opt.trueIota= obj.trueIota;
            opt.epsilon= obj.epsilon;
            opt.mask=obj.mask;

            maskIdx = find(obj.mask~=0);
            wvltIdx = find(obj.maskk~=0);
            opt.trueAlpha=obj.trueImg(maskIdx);
            
            %Sampling operator
            W=@(z) midwt(z,obj.wav,obj.dwt_L);
            Wt=@(z) mdwt(z,obj.wav,obj.dwt_L);

            obj.Psi=@(s) PsiFunc(s,W,obj.imgSize,maskIdx,wvltIdx);
            obj.Psit=@(s) PsitFunc(s,Wt,obj.imgSize,maskIdx,wvltIdx);
            fprintf('Configuration Finished!\n');
        end

        function nufftOps(obj)
            m_2D=[obj.imgSize, obj.imgSize];
            J=[1,1]*3;                       % NUFFT interpolation neighborhood
            K=2.^ceil(log2(m_2D*2));         % oversampling rate
            theta = (0:obj.prjNum-1)*360/obj.prjFull;

            r=pi*linspace(-1,1-2/obj.prjWidth,obj.prjWidth)';
            xc=r*cos(theta(:)'*pi/180);
            yc=r*sin(theta(:)'*pi/180);
            om=[yc(:), xc(:)];
            st=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');
            st.Num_pixel=obj.prjWidth;
            st.Num_proj=obj.prjNum;

            maskIdx = find(obj.mask~=0);
            % Zero freq at f_coeff(prjWidth/2+1)
            f_coeff=designFilter('renliang1',obj.prjWidth,obj.Ts);
            obj.Phi=@(s) PhiFunc51(s,f_coeff,st,obj.imgSize,obj.Ts,maskIdx);
            obj.Phit=@(s) PhitFunc51(s,f_coeff,st,obj.imgSize,obj.Ts,maskIdx);
            obj.FBP=@(s) FBPFunc6(s,theta,obj.Ts);
        end
        function cpuFanParOps(obj)
            conf.n=obj.imgSize; conf.prjWidth=obj.prjWidth;
            conf.np=obj.prjNum; conf.prjFull=obj.prjFull;
            conf.dSize=obj.dSize; %(n-1)/(Num_pixel+1);
            conf.effectiveRate=obj.effectiveRate;
            conf.d=obj.dist;

            mPrj(0,conf,'config');
            mPrj(0,0,'showConf');
            maskIdx = find(obj.mask~=0);
            obj.Phi =@(s) mPrj(maskFunc(s,maskIdx,conf.n),0,'forward')*obj.Ts;
            obj.Phit=@(s) maskFunc(mPrj(s,0,'backward'),maskIdx)*obj.Ts;
            obj.FBP = @(s) FBPFunc8(s,conf,obj.Ts,maskIdx)*obj.Ts;
        end
        function genOperators(obj)
            switch lower(obj.PhiMode)
                case 'basic'
                    nufftOps(obj);
                case lower('cpuPrj')
                    % can be cpu or gpu, with both fan or parallel projections
                    cpuFanParOps(obj);
                case lower('parPrj')
                    conf.bw=1; conf.nc=obj.imgSize; conf.nr=obj.imgSize; conf.prjWidth=obj.prjWidth;
                    conf.theta = (0:obj.prjNum-1)*360/obj.prjFull;
                    maskIdx = find(obj.mask~=0);
                    obj.Phi = @(s) mParPrj(s,maskIdx-1,conf,'forward')*obj.Ts;
                    obj.Phit= @(s) mParPrj(s,maskIdx-1,conf,'backward')*obj.Ts;
                    obj.FBP = @(s) FBPFunc7(s,obj.prjFull,obj.prjNum,obj.Ts,maskIdx)*obj.Ts;
                case 'weighted'
                    % Fessler's weighted methods
                    weight=exp(-y);
                    %weight=exp(-ones(size(y)));
                    weight=sqrt(weight/max(weight(:)));
                    y=y.*weight;

                    Phi=@(s) PhiFunc51(s,f_coeff,st,mx,Ts,maskIdx).*weight(:);
                    Phit=@(s) PhitFunc51(s.*weight(:),f_coeff,st,mx,...
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

                    Phi=@(s) PhiFunc2(s,f_coeff,stFwd,Num_pixel,Ts,maskIdx);
                    Phit=@(s) PhitFunc2(s,f_coeff,stFwd,Num_pixel,Ts,maskIdx);
                    FBP=Phit;
                otherwise
                    fprintf('Wrong mode for PhiMode: %s\n',PhiMode);
                    return;
            end
        end
        function junk(obj)
            Mask=double(Mask~=0);
            %figure; showImg(Mask);
            wvltIdx=find(maskk~=0);
            p_I=length(wvltIdx);
            maskIdx=find(Mask~=0);
            p_M=length(maskIdx);

            fprintf('Generating Func handles...\n');
            m=imgSize^2;

            H=@(s) Phi(Psi(s));
            Ht=@(s) Psit(Phit(s));

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


        end

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

            % before this, make sure max(CTdata)==1
            temp=size(obj.CTdata,1);
            obj.CTdata=[zeros(ceil((obj.prjWidth-temp)/2),obj.prjNum);...
                obj.CTdata;...
                zeros(floor((obj.prjWidth-temp)/2),obj.prjNum)];
            obj.y=-log(obj.CTdata(:)/max(obj.CTdata(:)));
        end

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

