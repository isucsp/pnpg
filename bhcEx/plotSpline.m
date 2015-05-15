
clear;
opt.E=25; opt.spectBasis='b1'; opt.logspan=3;
[y,args] = genBeamHarden({'Fe'},[],[],'showImg',false);
iota = args.iota(:);
epsilon = args.epsilon(:);
kappa = args.kappa(:);

iota=iota/sum(iota);

trueKappa=logspace(log10(min(kappa)),log10(max(kappa)),opt.E);
opt.E=length(trueKappa);
q=trueKappa(2)/trueKappa(1);

polymodel=Spline(opt.spectBasis,[trueKappa(1)/q; trueKappa(:); trueKappa(end)*q]);
polyIout = polymodel.polyIout;
polymodel.setPlot(kappa,iota,epsilon);

[upkappa,upiota]=getUpiota(epsilon,kappa,iota);
trueIe=interp1(log10(upkappa),upiota,log10(trueKappa(:)),'spline');
% there will be some points interplated negative and need to be removed
trueIe=max(trueIe,0);

forSave=[trueKappa(:), trueIe(:)];
save('trueI.data','forSave','-ascii');
forSave=[upkappa(:), upiota(:)];
save('continuousI.data','forSave','-ascii');

figure;
semilogx(upkappa,upiota,'b-'); hold on;
semilogx(trueKappa,trueIe,'r*-'); hold off;
xlabel('mass attenuation coefficient \kappa');
ylabel('density');
title('B1-spline approximation of \iota(\kappa)');

