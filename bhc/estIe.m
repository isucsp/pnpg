function [Ie,out] = estIe(kappa,Imea,PhiAlpha,spectBasis,noiseType)
    if(nargin==0)
        [a,~]=dbstack('-completenames'); a=a(1);
        [pathstr,~,~]=fileparts(a.file);
        addpath([pathstr '/../bhcEx']);
        clear a pathstr

        opt.beamharden=true; opt.errorType=0; opt.spectBasis='b1';
        opt.maxIeSteps=1e3; opt.stepShrnk=0.5; opt.thresh=1e-10;
        opt.prjFull = 360; opt.prjNum = 360; opt.snr=1e4;
        opt.noiseType='Poisson';
        opt.E=100; opt.logspan=3;

        [y,Phi,Phit,Psi,Psit,opt,FBP]=loadYang(opt);
        initSig = maskFunc(FBP(y),opt.mask~=0);

        Imea=0;
        temp = [opt.epsilon(1);(opt.epsilon(1:end-1)+opt.epsilon(2:end))/2;opt.epsilon(end)];
        deltaEpsilon = temp(2:end)-temp(1:end-1);
        for i=1:length(opt.epsilon)
            PhiAlpha = Phi(opt.trueAlpha) * opt.kappa(i); 
            Imea=Imea+exp(-PhiAlpha)*opt.iota(i)*deltaEpsilon(i);
        end
        Imea = Imea/max(Imea(:));
        Imea=exp(-y);

        kappa=logspace(-floor(opt.E/2)/(opt.E-1)*opt.logspan,...
            floor(opt.E/2-0.5)/(opt.E-1)*opt.logspan,opt.E);
        opt.E=length(kappa);

        [opt.upkappa,opt.upiota]=getUpiota(opt.epsilon,opt.kappa,opt.iota);
        opt.trueIe=interp1(log10(opt.upkappa),opt.upiota,log10(kappa(:)),'spline');
        % there will be some points interplated negative and need to be removed
        opt.trueIe=max(opt.trueIe,0);

        PhiAlpha=Phi(opt.trueAlpha);
        PhiFbp=Phi(initSig);

        [Ie,out] = estIe(kappa,Imea,Phi(opt.trueAlpha),opt.spectBasis,opt.noiseType);
        out.opt=opt;

        IelenRange=length(kappa):-1:2;

        return;
    end

    currNumIe=length(kappa);
    Ie=ones(size(kappa));
    maxIeSteps=5e3;
    switch lower(noiseType)
        case lower('Gaussian')
            IeStepFunc = @(A,III) gaussLI(Imea,A,III);
            hessian = @(A,III) getHess(A,(A*III).^2);
        case lower('Poisson')
            IeStepFunc = @(A,III) poissLI(Imea,A,III);
            hessian = @(A,III) getHess(A,(A*III));
    end
    out.Ie=[];
    out.GML=[];
    IelenRange=currNumIe:-1:2;
    ii=0; strlen=0;
    for Ielen=IelenRange
        ii=ii+1;
        temp=logspace(log10(kappa(1)),log10(kappa(end)),Ielen);
        Ie=interp1(log(kappa),Ie,log(temp),'spline');

        kappa=temp(:);
        temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
        polymodel = Spline(spectBasis,temp);
        polyIout = polymodel.polyIout;

        A = polyIout(PhiAlpha,[]);
        IeStep = NPG_ind(Ie,true,[],[],maxIeSteps);
        IeStep.func = @(III) IeStepFunc(A,III);
        Ie=IeStep.main();

        out.Ie{ii}=Ie;
        out.kappa{ii}=kappa;
        out.fisher{ii}=hessian(A,Ie);
        out.GML(ii)=IeStep.cost+0.5*log(det(out.fisher{ii}(Ie~=0,Ie~=0)));
        out.likelihood(ii)=IeStep.cost;
        out.rankA(ii)=rank(A);

        str=sprintf('|Ie|=%d, L=%g, |I|=%g, GML=%g',Ielen,out.likelihood(ii),det(out.fisher{ii}(Ie~=0,Ie~=0)),out.GML(ii));
        fprintf([repmat('\b',1,strlen) '%s'],str);
        strlen=length(str);
    end
    fprintf('\n');
end

function h = getHess(A,Iout)
    nc=size(A,2);
    h=A;
    for i=1:nc
        h(:,i)=h(:,i)./Iout;
    end
    h=A'*h;
end

