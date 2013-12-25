
close all;
clear;
fs=22;

loadMac

baseMaterial=26;
temp=mac{baseMaterial};
epsilon=logspace(-3,log10(20),4096);
epsilon=epsilon(:);
delta=temp(2:end,1)-temp(1:end-1,1);
idx=find(delta==0);
temp(idx,1)=temp(idx,1)-eps;
baseMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
original = baseMu;
[~, idx]=sort(baseMu,'descend');
%epsilon=epsilon(idx);
%loglog(mac{baseMaterial}(:,1)*1000,mac{baseMaterial}(:,2),'*');
loglog(epsilon*1e3,baseMu);
set(gca,'fontsize',fs);
hold on;

temp=mac{29};
delta=temp(2:end,1)-temp(1:end-1,1);
idx1=find(delta==0);
temp(idx1,1)=temp(idx1,1)*0.99;
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
figure;
semilogx(baseMu(idx), newMu(idx),'.');
hold on;
loglog(baseMu, newMu,'r');

hold on;
loglog(epsilon*1e3, newMu, 'r');
hold off;
xlabel('$\varepsilon$ (keV)','interpreter','latex');
ylabel('$\mu(\varepsilon)$ (cm$^2$/g)','interpreter','latex','fontname','times');
h=legend('Fe','Cu');
grid on
set(h,'interpreter','latex','fontname','times');
forsave=[epsilon*1e3, baseMu, newMu];
save('macFeCu.data','forsave','-ascii');
saveas(gcf,'macFeCu.eps','psc2'); return;
newMu=newMu./baseMu;
newMu=newMu; %-min(newMu(:))+eps;

%loglog(epsilon,newMu,'k-');
%figure(21); loglog(baseMu(idx),newMu(idx),'-'); hold on;
%figure(22); plot(baseMu(idx),newMu(idx),'-'); hold on;
%figure(23); semilogy(baseMu(idx),newMu(idx),'-'); hold on;
%figure(24); semilogy(epsilon,newMu,'-'); hold on;
figure(26); semilogx(epsilon(epsilon>0.01)*1000,newMu(epsilon>0.01),'-');
set(gca,'fontsize',fs);
xlabel('$\varepsilon$ (keV)','interpreter','latex');
ylabel('$\mu_{\rm{bone}}/\mu_{\rm{iron}}$','interpreter','latex','fontname','times');
fname='muRatioFeCuVsEpsilon';
saveas(gcf,[fname '.eps'],'psc2');

baseMu=baseMu(idx); newMu=newMu(idx);
idx=find(baseMu<30);
figure(25); semilogx(baseMu(idx),newMu(idx),'-'); hold on;
set(gca,'fontsize',fs);
xlabel('$\mu_{\rm{iron}}$ (cm$^2$/g)','interpreter','latex');
ylabel('$\mu_{\rm{bone}}/\mu_{\rm{iron}}$','interpreter','latex','fontname','times');
saveas(gcf,'muRatioFeCu.eps','psc2');

