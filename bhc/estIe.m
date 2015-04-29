function [Ie,out] = estIe(kappa,Imea,PhiAlpha,spectBasis,noiseType)
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
    ii=0;
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
        out.GML(ii)=IeStep.cost+0.5*log(det(hessian(A,Ie)));
        out.likelihood(ii)=IeStep.cost;
        out.rankA(ii)=rank(A);
    end
end

function h = getHess(A,Iout)
    nc=size(A,2);
    h=A;
    for i=1:nc
        h(:,i)=h(:,i)./Iout;
    end
    h=A'*h;
end

