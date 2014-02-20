function [CTdata, args] = genBeamHarden(varargin)
    %   Simulate and generate measurements for beamhardening effect
    %   Parameters includs:
    %   'symbol'
    %   'trueImg'
    %   'spark'
    %   'epsilon'
    %   'iota'
    %   'showImg'
    %   'saveMat'
    %   'filename'
    %
    %   Author: Renliang Gu (renliang@iastate.edu)
    %   $Revision: 0.1 $ $Date: Wed 05 Feb 2014 03:09:33 PM CST

    args = parseInputs(varargin{:});

    Phi = args.operators.Phi;
    %Phit = args.operators.Phit;
    FBP = args.operators.FBP;
    epsilon = args.epsilon(:);

    loadXrayMassCoef
    for i=1:length(symbols)
        if(strcmpi(args.symbol,symbols{i}))
            massAttenCoef=mac{i};
            %density = zaid(i,3);
            break;
        end
    end
    massAttenCoef(:,1)=massAttenCoef(:,1)*1e3;
    % from the original to 130 points, increase 19%
    % from 130 points to 1300 points, increasing 1e-6, which is not necessary.
    idx = massAttenCoef(:,1)>=min(epsilon) & massAttenCoef(:,1)<=max(epsilon);
    idx1= find(idx, 1 )-1; idx2 = find(idx, 1, 'last' )+1;
    kappa=interp1(massAttenCoef(idx1:idx2,1),massAttenCoef(idx1:idx2,2),...
        epsilon,'spline');

    if(args.showImg)
        %figure; loglog(massAttenCoef(:,1),massAttenCoef(:,2));
        %figure; semilogy(epsilon,kappa,'-');
        %hold on;loglog(massAttenCoef(:,1),massAttenCoef(:,2),'*');
        figure; loglog(massAttenCoef(idx1:idx2,2),...
            gampdf((massAttenCoef(idx1:idx2,1)-20)*16/100,5,1),'*');
        hold on; semilogy(kappa,args.iota);
        %figure; semilogx(massAttenCoef(:,1),args.iota);
        figure; semilogy(epsilon,args.iota,'*-');
    end

    epsilon=epsilon(:);
    deltaEpsilon=mean([epsilon(:) [epsilon(2:end); epsilon(end)]],2)-...
        mean([epsilon(:) [epsilon(1); epsilon(1:end-1)]],2);
    PhiAlpha=Phi(args.trueImg);

    Imea=0;
    for i=1:length(epsilon)
        Imea=Imea+exp(-PhiAlpha*kappa(i))*args.iota(i)*deltaEpsilon(i);
    end

    temp = max(Imea(:));
    Imea = Imea/temp;
    args.iota = args.iota/temp;

    args.kappa = kappa;
    args.s = PhiAlpha;
    CTdata=reshape(Imea,1024,[]);
    if(args.saveMat) save(filename,'CTdata'); end

    if(args.showImg)
        y=-log(Imea);
        rec=FBP(y);
        rec2=FBP(Phi(args.trueImg));
        figure; showImg(rec); figure; showImg(rec2);
    end
end

function args = parseInputs(varargin)
    args.symbol = 'Fe';
    args.spark = true;
    args.epsilon = (20:1:150)'; %keV
    args.showImg = true;
    args.saveMat = false;
    args.prjNum = 180;      % number of projections from degree at 0
    args.prjFull = 360;
    PhiMode = 'cpuPar'; %'nufft'; %'gpu'; %'cpu'

    for i=1:2:length(varargin)
        eval(['args.' varargin{i} ' = varargin{i+1};']);
    end

    if(~isfield(args,'trueImg'))
        args.trueImg=double(imread('binaryCasting.bmp'));
    end
    if(~isfield(args,'operators'))
        Ts = 0.008;
        theta = (0:args.prjNum-1)'/args.prjFull*360; % theta's in degree
        switch lower(args.PhiMode)
            case 'basic'
                n=1024;         % size of image
                prjWidth = 1024;
                m_2D=[n, n];
                J=[1,1]*3;                       % NUFFT interpolation neighborhood
                K=2.^ceil(log2(m_2D*2));         % oversampling rate

                r=pi*linspace(-1,1-2/prjWidth,prjWidth)';
                xc=r*cos(theta'*pi/180);
                yc=r*sin(theta'*pi/180);
                om=[yc(:), xc(:)];
                st=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');
                st.Num_pixel=prjWidth;
                st.Num_proj=args.prjNum;

                % Zero freq at f_coeff(prjWidth/2+1)
                f_coeff=designFilter('renliang1',prjWidth,Ts);
                args.operators.Phi=@(s) PhiFunc51(s,f_coeff,st,n,Ts);
                args.operators.Phit=@(s) PhitFunc51(s,f_coeff,st,n,Ts);
                args.operators.FBP=@(s) FBPFunc6(s,theta,Ts);
            case lower('cpuPar')
                conf.bw=1; conf.nc=1024; conf.nr=1024; conf.prjWidth=1024;
                conf.theta=theta;
                maskIdx=1:numel(args.trueImg);
                args.operators.Phi =@(s) mParPrj(s,maskIdx-1,conf,'forward')*Ts;
                args.operators.Phit=@(s) mParPrj(s,maskIdx-1,conf,'backward')*Ts;
                args.operators.FBP =@(s) FBPFunc7(s,args.prjFull,args.prjNum,Ts,maskIdx)*Ts;
            case lower('cpuFanPar')
                conf.n=1024; conf.prjWidth=1024;
                conf.np=args.prjNum; conf.prjFull=args.prjFull;
                conf.dSize=1; %(n-1)/(Num_pixel+1);
                conf.effectiveRate=1;
                conf.d=0;

                mPrj(0,conf,'config');
                maskIdx=1:numel(args.trueImg);
                args.operators.Phi =@(s) mPrj(maskFunc(s,maskIdx,conf.n),0,'forward')*Ts;
                args.operators.Phit=@(s) maskFunc(mPrj(s,0,'backward'),maskIdx)*Ts;
                args.operators.FBP =@(s) FBPFunc7(s,args.prjFull,args.prjNum,Ts,maskIdx)*Ts;
                %args.operators.FBP =@(s) FBPFunc8(s,conf,obj.Ts,maskIdx)*obj.Ts;
            otherwise
                fprintf('Wrong mode for PhiMode: %s\n',PhiMode);
                return;
        end
    end

    if(~isfield(args,'iota'))
        args.iota=gampdf((args.epsilon-20)*16/100,5,1);
        if(args.spark)
            args.iota(45)=args.iota(45)*1000;
            args.iota(22)=args.iota(22)*1000;
        end
    end

    if(~isfield(args,'filename'))
        if(args.spark)
            args.filename = 'CTdata_220SimSpark.mat';
        else
            args.filename = 'CTdata_220Sim.mat';
        end
    end
end
