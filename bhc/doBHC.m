function doBHC(y,distance,filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polychromatic Sparse Image Reconstruction and Mass Attenuation Spectrum 
%            Estimation via B-Spline Basis Function Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:  Changed to class oriented for easy configuration

    y=y-min(y(:));

    [N,M]=size(y);

    if(~exist('opt','var') || ~isfield(opt,'prjFull')) opt.prjFull=M; end
    if(~isfield(opt,'prjNum')) opt.prjNum=opt.prjFull; end

    if(mod(opt.prjFull,4)~=0)
        error(sprintf('prjFull=%d, is not dividable by 4\n',opt.prjFull));
    end
    if(mod(M,opt.prjFull)~=0)
        error(sprintf('prjFull=%d, doesn''t divide M=%d\n',opt.prjFull,M));
    end
    theta = (0:(opt.prjNum-1))*M/opt.prjFull +1;
    y=y(:,theta);

    conf=ConfigCT();

    daub = 2; dwt_L=6;        %levels of wavelet transform
    maskType='CircleMask';

    conf.PhiMode = 'gpuPrj'; %'parPrj'; %'basic'; %'gpuPrj'; %
    conf.imgSize = 2^floor(log2(N));
    conf.prjWidth = N;
    conf.prjFull = opt.prjFull;
    conf.prjNum = opt.prjNum;
    conf.dSize = conf.imgSize/conf.prjWidth;
    conf.effectiveRate = 1;
    conf.dist = distance/(conf.prjWidth)*(conf.imgSize);
    conf.Ts =1e-2;

    if(strcmpi(maskType,'CircleMask'))
        % reconstruction mask (which pixels do we estimate?)
        mask = Utils.getCircularMask(conf.imgSize);
        wvltName = sprintf('MaskWvlt%dCircleL%dD%d.mat',conf.imgSize,dwt_L,daub);
        if(exist(wvltName,'file'))
            load(wvltName);
        else
            maskk=wvltMask(mask,dwt_L,daub,wvltName);
        end
    elseif(strcmpi(maskType,'none'))
        mask = []; maskk=[];
    end
    [Phi,Phit,FBP]=conf.genOperators(mask);
    [Psi,Psit]=Utils.getPsiPsit(daub,dwt_L,mask,maskk);
    opt.mask=mask; opt.maskk=maskk;

    fprintf('Configuration Finished!\n');

    opt.beamharden=true; opt.spectBasis='dis'; opt.E=20;
    opt.maxItr=2; opt.thresh=1e-6;

    initSig = maskFunc(FBP(y),opt.mask~=0);

    i=1; j=1;
    fprintf('%s, i=%d, j=%d\n','Filtered Backprojection',i,j);
    fbp.img=FBP(y);
    fbp.alpha=fbp.img(opt.mask~=0);

    Oopt=opt;

    % unknown ι(κ), NPG-AS
    for j=5:-1:2
        fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
        u  =  10^(-5);
        opt=Oopt; opt.u=u*10^(j-3); opt.proximal='tvl1';
        if(j==5)
            npgTV{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
        else
            opt.Ie=npgTV{i,j+1}.Ie;
            npgTV{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,...
                npgTV{i,j+1}.alpha,opt);
        end

        opt=Oopt; opt.u=u*10^(j-3); opt.proximal='wvltADMM';
        npgWV{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
        if(j==5)
            npgWV{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
        else
            opt.Ie=npgWV{i,j+1}.Ie;
            npgWV{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,...
                npgWV{i,j+1}.alpha,opt);
        end

        save(filename);
    end
end

