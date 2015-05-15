
clear

opt.beamharden=true; opt.errorType=0; opt.spectBasis='b1';
opt.noiseType='Poisson'; opt.prjFull=360; opt.prjNum=360;
opt.E=35; opt.logspan=3;

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
q=kappa(2)/kappa(1);
polymodel=Spline(opt.spectBasis,[kappa(1)/q; kappa(:); kappa(end)*q]);
polymodel.setPlot(opt.kappa,opt.iota,opt.epsilon);

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

switch lower(opt.noiseType)
    case lower('Gaussian')
        IeStepFunc = @(AA,III) gaussLI(Imea,AA,III);
    case lower('Poisson')
        IeStepFunc = @(AA,III) poissLI(Imea,AA,III);
end

fprintf('cost=%e\n',IeStepFunc(trueA,opt.trueIe));

s=linspace(min(PhiAlpha),max(PhiAlpha),100);
idx=randi(length(PhiAlpha),1000,1);

figure;
plot(PhiAlpha(idx),-log(Imea(idx)),'.'); hold on;
plot(PhiFbp(idx),-log(Imea(idx)),'g.');
plot(s,-log(polyIout(s,opt.trueIe)),'r-');
legend('perfect reconstruction', 'FBP reconstruction', 'fitted curve');
xlabel('\Phi\alpha');
ylabel('I^{out}=-ln[ \int \iota(\kappa) exp( -\kappa\Phi\alpha ) d\kappa  ]');

forSave=[PhiAlpha(idx),-log(Imea(idx))];
save('test1.data','forSave','-ascii');
forSave=[PhiFbp(idx),-log(Imea(idx))];
save('test2.data','forSave','-ascii');
forSave=[s(:), -log(polyIout(s,opt.trueIe))];
save('test3.data','forSave','-ascii');

!for i in `seq 1 3`; do echo "" >> test$i.data; done
!cat test[1-3].data > linearization.data
!rm test[1-3].data

