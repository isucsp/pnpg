function [Ie,out] = testIeStep()
    opt.beamharden=true; opt.errorType=0; opt.spectBasis='dis';
    opt.maxIeSteps=1e4; opt.stepShrnk=0.5; opt.thresh=1e-10;
    opt.prjFull = 360; opt.prjNum = 360; opt.snr=1e4;
    opt.noiseType='Poisson';
    opt.E=50; opt.logspan=3;

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

    [Ie,kappa,polymodel,opt] = setPolyModel(opt);
    polyIout = polymodel.polyIout;

    [opt.upkappa,opt.upiota]=getUpiota(opt.epsilon,opt.kappa,opt.iota);
    opt.trueIe=interp1(log10(opt.upkappa),opt.upiota,log10(kappa(:)),'spline');
    % there will be some points interplated negative and need to be removed
    opt.trueIe=max(opt.trueIe,0);
    trueA=polyIout(Phi(opt.trueAlpha),[]);
    norm(y+log(trueA*opt.trueIe))

    switch lower(opt.noiseType)
        case lower('Gaussian')
            IeStepFunc = @(A,III) gaussLI(Imea,A,III);
        case lower('Poisson')
            IeStepFunc = @(A,III) poissLI(Imea,A,III);
    end

    fprintf('cost=%e\n',IeStepFunc(trueA,opt.trueIe));

    B=eye(opt.E); b=zeros(opt.E,1);

    %Ie=opt.trueIe;
    IeStep=NPG_ind(Ie,true,polyIout(0,[]),1,1,opt.stepShrnk,opt.thresh);
    IeStep=NPG_ind(Ie,true,[],1,1,opt.stepShrnk,opt.thresh);
    IeStep.func=@(III) IeStepFunc(trueA,III);

    tic;
    strlen=0;

    for i=1:opt.maxIeSteps;
        IeStep.main();
        out.cost(i)=IeStep.cost;
        out.difIe(i)=IeStep.difIe;
        out.stepSize(i)=1/IeStep.t;
        if(i>1)
            out.difCost(i)=abs((out.cost(i)-out.cost(i-1)))/out.cost(i);
        end
        out.time(i)=toc;

        str=sprintf('i=%d, time=%g, cost=%e',i,out.time(i),IeStep.cost);
        fprintf([repmat('\b',1,strlen) '%s'],str);
        strlen=length(str);
    end
    Ie = IeStep.Ie;
    out.kappa=kappa;
    out.opt=opt;
    fprintf('\n');
end

function [Ie,kappa,polymodel,opt] = setPolyModel(opt)
    temp=logspace(-floor(opt.E/2)/(opt.E-1)*opt.logspan,...
        floor(opt.E/2-0.5)/(opt.E-1)*opt.logspan,opt.E);
    opt.E=length(temp);
    Ie=zeros(opt.E,1); Ie(floor(opt.E/2)+1)=1;

    kappa=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(opt.K-1)*eps);
    temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
    polymodel = Spline(opt.spectBasis,temp);
    polymodel.setPlot(opt.kappa,opt.iota,opt.epsilon);

    Ie=Ie/polymodel.polyIout(0,Ie);
end

