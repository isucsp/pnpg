
clear;

loadMac

%baseMaterial=bone.mac;
temp=bone.mac; %mac{baseMaterial};
temp1=temp;
epsilon=logspace(-3,log10(20),4096);
epsilon=epsilon(:);
delta=temp(2:end,1)-temp(1:end-1,1);
idx=find(delta==0);
temp(idx,1)=temp(idx,1)*0.999;
baseMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
original = baseMu;
[~, idx]=sort(baseMu,'descend');
%epsilon=epsilon(idx);
loglog(temp1(:,1),temp1(:,2),'*');
hold on;
loglog(epsilon,baseMu);
hold off;

temp=brain.mac;
delta=temp(2:end,1)-temp(1:end-1,1);
idx1=find(delta==0);
temp(idx1,1)=temp(idx1,1)*0.999;
newMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
figure;
loglog(temp(:,1),temp(:,2),'*'); hold on;
loglog(epsilon,newMu);
hold off;
%figure; plot(baseMu, newMu); hold on; plot(baseMu,baseMu,'r');
figure; loglog(baseMu, newMu); hold on; plot(baseMu,baseMu,'r');
newMu=newMu./baseMu; newMu=newMu; %-min(newMu(:))+eps;
figure;
loglog(baseMu(idx), newMu(idx),'.');
hold on;
loglog(baseMu, newMu,'r');

fs=30;
%loglog(epsilon,newMu,'k-');
%figure(21); loglog(baseMu(idx),newMu(idx),'-'); hold on;
%figure(22); plot(baseMu(idx),newMu(idx),'-'); hold on;
%figure(23); semilogy(baseMu(idx),newMu(idx),'-'); hold on;
%figure(24); semilogy(epsilon,newMu,'-'); hold on;
figure(26); semilogx(epsilon(epsilon>0.01)*1000,newMu(epsilon>0.01),'-');
set(gca,'fontsize',fs);
xlabel('$\varepsilon$ (keV)','interpreter','latex');
ylabel('$\kappa_{\rm{bone}}/\kappa_{\rm{iron}}$','interpreter','latex','fontname','times');
fname='kappaRatioVSEpsilon';
saveas(gcf,[fname '.eps'],'psc2');

%baseMu=baseMu(idx); newMu=newMu(idx);
idx=find(baseMu<50);
figure(25); semilogx(baseMu(idx),newMu(idx),'-'); hold on;
temp=[baseMu(idx), newMu(idx)];
save('kappaRatioBrainVsBone.data','temp','-ascii');
set(gca,'fontsize',fs-5);
xlabel('$\kappa_{\rm{iron}}$ (cm$^2$/g)','interpreter','latex','fontsize',fs);
ylabel('$\kappa_{\rm{bone}}/\kappa_{\rm{iron}}$','interpreter','latex','fontname','times','fontsize',fs);
saveas(gcf,'kappaRatioVSIron.eps','psc2');
close all

