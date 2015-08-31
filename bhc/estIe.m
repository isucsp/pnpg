function [Ie,out] = estIe(Erange,Imea,PhiAlpha,spectBasis,noiseType,upkappa,upiota)
if(nargin==0)
    [a,~]=dbstack('-completenames'); a=a(1);
    [pathstr,~,~]=fileparts(a.file);
    addpath([pathstr '/../bhcEx']);
    clear a pathstr

    opt.beamharden=true; opt.errorType=0; opt.spectBasis='b1';
    opt.maxIeSteps=5e3; opt.stepShrnk=0.5; opt.thresh=1e-10;
    opt.prjFull = 360; opt.prjNum = 360; opt.snr=1e4;
    opt.logspan=3;

    opt.noiseType='';
    [y,Phi,Phit,Psi,Psit,opt,FBP]=loadYang(opt);
    Imea=exp(-y);
    PhiAlpha=Phi(opt.trueAlpha);
    [opt.upkappa,opt.upiota]=getUpiota(opt.epsilon,opt.kappa,opt.iota);

    opt.noiseType='Poisson';
    [Ie,out]=estIe(17,Imea,PhiAlpha,opt.spectBasis,opt.noiseType,opt.upkappa,opt.upiota);
    out.opt=opt;

    return;
end

    logspan=3;
    Ielen=min(Erange);
    kappa=logspace(-floor(Ielen/2)/(Ielen-1)*logspan, floor(Ielen/2-0.5)/(Ielen-1)*logspan,Ielen)';
    %Ie=ones(Ielen,1);
    Ie=zeros(Ielen,1); Ie(floor(Ielen/2)+1)=1;
    maxIeSteps=5e3;
    switch lower(noiseType)
        case lower('Gaussian')
            IeStepFunc = @(A,III) gaussLI(Imea,A,III);
            hessian = @(A,III) getHess(A,(A*III).^2);
        case lower('Poisson')
            IeStepFunc = @(A,III) poissLI(Imea,A,III);
            hessian = @(A,III) getHess(A,(A*III));
    end
    IeStep=LBFGSB(Ie,1e3);
    IeStep.thresh=1e-19;
    IeStep.opts.printEvery=100;

    out.Ie=[];
    out.GML=[];
    ii=0; strlen=0;
    for Ielen=Erange
        ii=ii+1;

        temp=logspace(-floor(Ielen/2)/(Ielen-1)*logspan, floor(Ielen/2-0.5)/(Ielen-1)*logspan,Ielen)';
        Ie=interp1(log(kappa),Ie,log(temp),'spline'); Ie(Ie<0)=0;
        kappa=temp;

        temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
        polymodel = Spline(spectBasis,temp);
        polyIout = polymodel.polyIout;

        A = polyIout(PhiAlpha,[]);
        IeStep.Ie=Ie/polyIout(0,Ie);
        IeStep.func = @(III) IeStepFunc(A,III);
            IeStep.opts.maxIts
        Ie=IeStep.main();

        out.Ie{ii}=Ie;
        out.kappa{ii}=kappa;
        out.fisher{ii}=hessian(A,Ie);
        out.GML(ii)=IeStep.cost+0.5*log(det(out.fisher{ii}(Ie~=0,Ie~=0)));
        out.NLL(ii)=IeStep.cost;
        out.rankA(ii)=rank(A);

        str=sprintf('|Ie|=%d, L=%g, |I|=%g, GML=%g',Ielen,out.NLL(ii),det(out.fisher{ii}(Ie~=0,Ie~=0)),out.GML(ii));
        fprintf([repmat('\b',1,strlen) '%s\n'],str);
        strlen=length(str)*0;
        figure; semilogx(kappa,Ie,'*-');
        if exist('upkappa','var')
            hold on; semilogx(upkappa,upiota,'r'); hold off;
        end
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

