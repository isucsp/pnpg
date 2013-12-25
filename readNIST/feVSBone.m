
clear;

loadMac

baseMaterial=26;
temp=mac{baseMaterial};
epsilon=logspace(-3,log10(20),4096);
epsilon=epsilon(:);
delta=temp(2:end,1)-temp(1:end-1,1);
idx=find(delta==0);
temp(idx,1)=temp(idx,1)*0.999;
baseMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
original = baseMu;
[~, idx]=sort(baseMu,'descend');
%epsilon=epsilon(idx);
loglog(mac{baseMaterial}(:,1),mac{baseMaterial}(:,2),'*');
hold on;
loglog(epsilon,baseMu);
hold off;

temp=bone.mac; %mac{i};
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

fs=22;
%loglog(epsilon,newMu,'k-');
%figure(21); loglog(baseMu(idx),newMu(idx),'-'); hold on;
%figure(22); plot(baseMu(idx),newMu(idx),'-'); hold on;
%figure(23); semilogy(baseMu(idx),newMu(idx),'-'); hold on;
%figure(24); semilogy(epsilon,newMu,'-'); hold on;
figure(26); semilogx(epsilon(epsilon>0.01)*1000,newMu(epsilon>0.01),'-');
set(gca,'fontsize',fs);
xlabel('$\varepsilon$ (keV)','interpreter','latex');
ylabel('$\mu_{\rm{bone}}/\mu_{\rm{iron}}$','interpreter','latex','fontname','times');
fname='muRatioVSEpsilon';
saveas(gcf,[fname '.eps'],'psc2');

baseMu=baseMu(idx); newMu=newMu(idx);
idx=find(baseMu<30);
figure(25); semilogx(baseMu(idx),newMu(idx),'-'); hold on;
set(gca,'fontsize',fs);
xlabel('$\mu_{\rm{iron}}$ (cm$^2$/g)','interpreter','latex');
ylabel('$\mu_{\rm{bone}}/\mu_{\rm{iron}}$','interpreter','latex','fontname','times');
saveas(gcf,'muRatioVSIron.eps','psc2');
