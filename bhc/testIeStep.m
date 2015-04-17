function testIeStep()
    opt.beamharden=true; opt.errorType=0; opt.spectBasis='dis';
    opt.maxIeSteps=1000; opt.stepShrnk=0.5;
    opt.prjFull = 360; opt.prjNum = 360; opt.snr=1e4;
    opt.E=17; opt.logspan=3;
    [y,Phi,Phit,Psi,Psit,opt,FBP]=loadYang(opt);
    initSig = maskFunc(FBP(y),opt.mask~=0);

    Imea=0;
    temp = [opt.epsilon(1);(opt.epsilon(1:end-1)+opt.epsilon(2:end))/2;opt.epsilon(end)];
    deltaEpsilon = temp(2:end)-temp(1:end-1);
    for i=1:length(opt.epsilon)
        PhiAlpha = Phi(opt.trueAlpha) * opt.trueKappa(i); 
        Imea=Imea+exp(-PhiAlpha)*opt.trueIota(i)*deltaEpsilon(i);
    end
    Imea = Imea/max(Imea(:));

    Imea=exp(-y);

    [Ie,kappa,polymodel,opt] = setPolyModel(opt);
    polyIout = polymodel.polyIout;

    temp1 = [opt.epsilon(1);(opt.epsilon(1:end-1)+opt.epsilon(2:end))/2;opt.epsilon(end)];
    temp2 = [opt.trueKappa(1);(opt.trueKappa(1:end-1)+opt.trueKappa(2:end))/2;opt.trueKappa(end)];
    temp1 = temp1(2:end)-temp1(1:end-1);
    temp2 = temp2(2:end)-temp2(1:end-1);
    opt.trueUpiota=abs(opt.trueIota.*temp1./temp2);

    opt.trueIe=interp1(opt.trueKappa, opt.trueUpiota,kappa(:),'spline');
    % there will be some points interplated negative and need to be removed
    opt.trueIe=max(opt.trueIe,0);
    trueA=polyIout(Phi(opt.trueAlpha),[]);
    norm(y+log(trueA*opt.trueIe))
    fprintf('cost=%e\n',gaussLI(Imea,trueA,opt.trueIe));

    trueA=polyIout(Phi(opt.trueAlpha*kappa(2)/kappa(1)),[]);
    norm(y+log(trueA*opt.trueIe))
    fprintf('cost=%e\n',gaussLI(Imea,trueA,opt.trueIe));

    B=eye(opt.E); b=zeros(opt.E,1);
    temp = polyIout(0,[]); B=[B; -temp(:)'/norm(temp)]; b=[b; -1/norm(temp)];

    IeStep = ActiveSet(B,b,Ie,opt.maxIeSteps,opt.stepShrnk);
    IeStep.func = @(III) gaussLI(Imea,trueA,III);
    IeStep.main();
    fprintf('cost=%e\n',IeStep.cost);

    IeStep = NPG_ind(Ie,true,[],1,opt.maxIeSteps,opt.stepShrnk);
    IeStep.func = @(III) gaussLI(Imea,trueA,III);
    IeStep.main();
    fprintf('cost=%e\n',IeStep.cost);
    
    keyboard;

    IeStep = NPG_ind(Ie,false,[],1,opt.maxIeSteps,opt.stepShrnk);
    IeStep.func = @(III) gaussLI(Imea,trueA,III);
    IeStep.main();
    fprintf('cost=%e\n',IeStep.cost);

    IeStep = NPG_ind(Ie,false,polyIout(0,[]),1,opt.maxIeSteps,opt.stepShrnk);
    IeStep.func = @(III) gaussLI(Imea,trueA,III);
    IeStep.main();
    fprintf('cost=%e\n',IeStep.cost);

    IeStep = NPG_ind(Ie,true,polyIout(0,[]),1,opt.maxIeSteps,opt.stepShrnk);
    IeStep.func = @(III) gaussLI(Imea,trueA,III);
    IeStep.main();
    fprintf('cost=%e\n',IeStep.cost);
end

function [Ie,kappa,polymodel,opt] = setPolyModel(opt)
    temp=logspace(-floor(opt.E/2)/(opt.E-1)*opt.logspan,...
        floor(opt.E/2-0.5)/(opt.E-1)*opt.logspan,opt.E);
    opt.E=length(temp);
    Ie=zeros(opt.E,1); Ie(floor(opt.E/2)+1)=1;

    kappa=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(opt.K-1)*eps);
    temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
    polymodel = Spline(opt.spectBasis,temp);
    polymodel.setPlot(opt.trueKappa,opt.trueIota,opt.epsilon);

    Ie=Ie/polymodel.polyIout(0,Ie);

end


