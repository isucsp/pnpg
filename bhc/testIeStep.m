function [Ie,out] = testIeStep()

    [a,~]=dbstack('-completenames'); a=a(1);
    [pathstr,~,~]=fileparts(a.file);
    addpath([pathstr '/../bhcEx']);
    clear a pathstr

    opt.beamharden=true; opt.errorType=0; opt.spectBasis='b1';
    opt.maxIeSteps=1e5; opt.stepShrnk=0.5; opt.thresh=1e-19;
    opt.prjFull = 360; opt.prjNum = 360; opt.snr=1e4;
    opt.noiseType='Poisson';
    opt.E=50; opt.logspan=2.5;

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

    PhiAlpha=Phi(opt.trueAlpha);
    PhiFbp=Phi(initSig);
    trueA=polyIout(PhiAlpha,[]);
    A=polyIout(PhiFbp,[]);
    disp([norm(y+log(trueA*opt.trueIe)) norm(y+log(A*opt.trueIe))]);

    keyboard

    switch lower(opt.noiseType)
        case lower('Gaussian')
            IeStepFunc = @(AA,III) gaussLI(Imea,AA,III);
        case lower('Poisson')
            IeStepFunc = @(AA,III) poissLI(Imea,AA,III);
    end

    fprintf('cost=%e\n',IeStepFunc(trueA,opt.trueIe));

    B=eye(opt.E); b=zeros(opt.E,1);

    Ie=Ie*0+1; Ie=Ie/polyIout(0,Ie);
    %Ie=opt.trueIe;
    IeStep=NPG_ind(Ie,true,polyIout(0,[]),1,1,opt.stepShrnk,opt.thresh);
    IeStep=NPG_ind(Ie,true,[],1,1,opt.stepShrnk,opt.thresh);
    IeStep.func=@(III) IeStepFunc(trueA,III);

    tic;
    strlen=0;

    s=linspace(min(PhiAlpha),max(PhiAlpha),1000);

    h1=figure; h2=figure;

    for i=1:opt.maxIeSteps;

        if(exist('h1','var'))
            idx=randi(length(PhiAlpha),1000,1);
            figure(h1); plot(PhiAlpha(idx),-log(Imea(idx)),'.'); hold on;
            plot(PhiFbp(idx),-log(Imea(idx)),'g.');
            plot(s,-log(polyIout(s,IeStep.Ie)),'r-'); hold off;
            figure(h2); semilogx(opt.upkappa,opt.upiota,'r'); hold on;
            plot(kappa,IeStep.Ie,'.-'); hold off;
            drawnow;
        end

        IeStep.main();
        out.cost(i)=IeStep.cost;
        out.difIe(i)=IeStep.difIe;
        out.stepSize(i)=1/IeStep.t;
        if(i>1)
            out.difCost(i)=abs((out.cost(i)-out.cost(i-1)))/out.cost(i);
        end
        out.time(i)=toc;

        str=sprintf('i=%d, time=%.2g, cost=%.10g, difIe=%g, stepSize=%g',i,out.time(i),IeStep.cost,IeStep.difIe,out.stepSize(i));
        if(i>1)
            str=sprintf('%s, difCost=%g',str,out.difCost(i));
        end
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

% Conclusions form this test file:
% - it is not good to initilize Ie by a pulse with its position centered
% - better to begin with eps*ones
