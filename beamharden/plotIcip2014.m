
style1 = {'b-','r-','b:','r:','b--','r--','-','*','.','+','d'};
style2 = {'b-*','r-*','b:o','r:o','b--s','r--s'};
filename='runNDE2014_1.mat';

%% NDE report plot
if(1)
    aa=0;
    legStr=[];
    toPlot=[];
    prj = [60, 80, 100, 120, 180, 360]/2;
    toPlot=[toPlot prj(:)];

    % unkown Ie, CPLS
    load(filename,'out012');
    out=out012;
    rse=[];
    for i=1:6; rse(i)=out{i,1}.RMSE(end); end
    toPlot=[toPlot rse(:)];
    temp=[];
    for i=1:6; temp(i)=out{i,4}.RMSE(end); end
    clear('out012');

    % known-Ie CPLS
    load(filename,'out001');
    figure;
    out=out001;
    rse=[out{1,5}.RMSE(end);
         out{2,4}.RMSE(end);
         out{3,4}.RMSE(end);
         out{4,4}.RMSE(end);
         out{5,4}.RMSE(end);
         out{6,4}.RMSE(end);];
    toPlot=[toPlot rse(:)];
    clear('out001');

    load(filename,'out003');
    out=out003;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=out{i}.RMSE;
    end
    toPlot=[toPlot rse(:)];
    clear('out003');

    load(filename,'out005');
    out=out005;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=out{i}.RMSE;
    end
    toPlot=[toPlot rse(:)];
    clear('out005');

    load(filename,'out002');
    out=out002;
    rse = ones(size(out))*inf;
    for i=1:size(out,1)
        for j=1:size(out,2)
            if(~isempty(out{i,j}))
                rse(i,j)=out{i,j}.RMSE(end);
            end
        end
    end
    toPlot=[toPlot min(rse,[],2)];
    clear('out002');

    load(filename,'out004');
    out=out004;
    rse = ones(size(out))*inf;
    for i=1:size(out,1)
        for j=1:size(out,2)
            if(~isempty(out{i,j}))
                rse(i,j)=out{i,j}.RMSE(end);
            end
        end
    end
    toPlot=[toPlot min(rse,[],2)];
    clear('out004');

    toPlot=[toPlot, temp(:)];

    save('RSEvsPrjNum.data','toPlot','-ascii');
end
return;

%% Compare the convergence speed between Old algorithm and the new for 180 parallel beam 
filename='runNDE2014_1.mat';
t=1; cost(:,t)=1:2000; rmse(:,t)=1:2000; time(:,t)=1:2000;

% FISTA, no restart, skipIe, no continuation,
load(filename,'out001');
out=out001{6,4}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;

load(filename,'out012');
% FISTA, restart, no skipIe, no continuation
out=out012{6,1}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;
% FISTA, restart, no skipIe, continuation
out=out012{6,2}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;
% FISTA, restart, skipIe, continuation
out=out012{6,3}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;
% FISTA, restart, skipIe, no continuation
out=out012{6,4}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;

load(filename,'out021');
% FISTA, restart, no skipIe, no continuation
out=out021{1,1}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;
% NCG, no skipIe, no continuation
out=out021{1,2}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;
% FISTA, restart, no skipIe, continuation
out=out021{2,1}; t=t+1;
cost(1:length(out.cost),t)=out.cost; rmse(1:length(out.RMSE),t)=out.RMSE; time(1:length(out.time),t)=out.time;

mincost=reshape(cost(:,2:end),[],1);
mincost=min(mincost(mincost>0));
% compare restart for FISTA, restart gives smaller cost, although rmse is
% not the minimum
figure; subplot(2,1,1); semilogy(cost(:,2)-mincost); hold on; plot(cost(:,6)-mincost,'r');
legend('no restart','restart'); title(sprintf('%g, %g',min(cost(cost(:,2)>0,2)),min(cost(cost(:,6)>0,6))));
subplot(2,1,2); semilogy(rmse(:,2)); hold on; plot(rmse(:,6),'r');
legend('no restart','restart'); title(sprintf('%g, %g',min(rmse(rmse(:,2)>0,2)),min(rmse(rmse(:,6)>0,6))));

% compare skipIe(restart, no continuation), skip to achieve better rmse,
% while keeping higher cost
figure; subplot(2,1,1); semilogy(cost(:,3)-mincost); hold on; plot(cost(:,6)-mincost,'r');plot(cost(:,7)-mincost,'g');legend('noskip','skip','noskip');
title(sprintf('%g, %g, %g',min(cost(cost(:,3)>0,3)),min(cost(cost(:,6)>0,6)),min(cost(cost(:,7)>0,7))));
subplot(2,1,2); semilogy(rmse(:,3)); hold on; plot(rmse(:,6),'r');plot(rmse(:,7),'g');legend('noskip','skip','noskip');
title(sprintf('%g, %g, %g',min(rmse(rmse(:,3)>0,3)),min(rmse(rmse(:,6)>0,6)),min(rmse(rmse(:,7)>0,7))));

% compare skipIe(restart, continuation), see above for conclusion
figure; subplot(2,1,1); semilogy(cost(:,4)-mincost); hold on; plot(cost(:,5)-mincost,'r');plot(cost(:,9)-mincost,'g');legend('noskip','skip','noskip');
title(sprintf('%g, %g, %g',min(cost(cost(:,4)>0,4)),min(cost(cost(:,5)>0,5)),min(cost(cost(:,9)>0,9))));
subplot(2,1,2); semilogy(rmse(:,4)); hold on; plot(rmse(:,5),'r');plot(rmse(:,9),'g');legend('noskip','skip','noskip');
title(sprintf('%g, %g, %g',min(rmse(rmse(:,4)>0,4)),min(rmse(rmse(:,5)>0,5)),min(rmse(rmse(:,9)>0,9))));

% compare continuation (restart, skipIe)
figure; subplot(2,1,1); semilogy(cost(:,5)-mincost); hold on; plot(cost(:,6)-mincost,'r');
legend('conti','no conti'); title(sprintf('%g, %g',min(cost(cost(:,5)>0,5)),min(cost(cost(:,6)>0,6))));
subplot(2,1,2); semilogy(rmse(:,5)); hold on; plot(rmse(:,6),'r');
legend('conti','no conti'); title(sprintf('%g, %g',min(rmse(rmse(:,5)>0,5)),min(rmse(rmse(:,6)>0,6))));

% compare continuation(restart, no skipIe)
figure; subplot(2,1,1); semilogy(cost(:,4)-mincost); hold on; plot(cost(:,9)-mincost,'g');plot(cost(:,3)-mincost,'r');plot(cost(:,7)-mincost,'k');legend('conti','conti','no conti','no conti');
title(sprintf('%g, %g, %g, %g',min(cost(cost(:,4)>0,4)),min(cost(cost(:,9)>0,9)),min(cost(cost(:,3)>0,3)),min(cost(cost(:,7)>0,7))));
subplot(2,1,2); semilogy(rmse(:,4)); hold on; plot(rmse(:,9),'g');plot(rmse(:,3),'r');plot(rmse(:,7),'r');legend('conti','conti','no conti','no conti');
title(sprintf('%g, %g, %g, %g',min(rmse(rmse(:,4)>0,4)),min(rmse(rmse(:,9)>0,9)),min(rmse(rmse(:,3)>0,3)),min(rmse(rmse(:,7)>0,7))));

% Compare the methods
figure; subplot(2,1,1); semilogy(cost(:,7)-mincost); hold on; plot(cost(:,8)-mincost,'r');
legend('conti','no conti'); title(sprintf('%g, %g',min(cost(cost(:,7)>0,7)),min(cost(cost(:,8)>0,8))));
subplot(2,1,2); semilogy(rmse(:,7)); hold on; plot(rmse(:,8),'r');
legend('conti','no conti'); title(sprintf('%g, %g',min(rmse(rmse(:,7)>0,7)),min(rmse(rmse(:,8)>0,8))));
figure; subplot(2,1,1); semilogy(time(:,7),cost(:,7)-mincost); hold on; plot(time(:,8),cost(:,8)-mincost,'r');
legend('conti','no conti'); title(sprintf('%g, %g',min(cost(cost(:,7)>0,7)),min(cost(cost(:,8)>0,8))));
subplot(2,1,2); semilogy(time(:,7),rmse(:,7)); hold on; plot(time(:,8),rmse(:,8),'r');
legend('conti','no conti'); title(sprintf('%g, %g',min(rmse(rmse(:,7)>0,7)),min(rmse(rmse(:,8)>0,8))));

%% plot the the estimated spectrum, images and everything for display
clear;
filename='runNDE2014_1.mat';
load(filename,'out012');
i=2;
out=out012{i,1};
polymodel = Spline(out.opt.spectBasis,out.kappa);
polymodel.setPlot(out.opt.trueKappa,out.opt.trueIota,out.opt.epsilon);
[tK,tU,kappa,Ie]=polymodel.plotSpectrum(out.Ie);
%loglog(tK,tU); hold on; loglog(kappa/1.81,Ie,'r*-'); hold off;
temp=[tK,tU];
save('continuousI.data','temp','-ascii');
temp=[kappa(Ie>0),Ie(Ie>0)];
save('estI.data','temp','-ascii');
[tK,tU,kappa,Ie]=polymodel.plotSpectrum(out012{i,4}.Ie);
temp=[kappa,Ie];
save('trueI.data','temp','-ascii');
fprintf('FISTA_ADMM: %g\n',out.RMSE(end));
img=showImgMask(out.alpha,out.opt.mask);
imwrite(img/max(img(:)),'FISTA_ADMM.png','png');
clear out012;

load(filename,'out002'); out=out002; o2=showResult(out,2,'RMSE');
j=find(o2(i,:)==min(o2(i,o2(i,:)>0)));
out=out002{i,j};
fprintf('FPCAS: %g\n',out.RMSE(end));
img=showImgMask(out.alpha,out.opt.mask);
imwrite(img/max(img(:)),'FPCAS.png','png');

load(filename,'out003'); out=out003{i};
fprintf('FBP: %g\n',out.RMSE(end));
img=showImgMask(out.alpha,out002{i,j}.opt.mask);
imwrite(img/max(img(:)),'FBP.png','png');
load(filename,'out005'); out=out005{i};
fprintf('FBPwLin: %g\n',out.RMSE(end));
img=showImgMask(out.alpha,out002{i,j}.opt.mask);
imwrite(img/max(img(:)),'FBPwLin.png','png');
clear out003, out005; clear out002;

load(filename,'out004'); out=out004; o2=showResult(out,2,'RMSE');
j=find(o2(i,:)==min(o2(i,o2(i,:)>0)));
out=out004{i,j};
fprintf('FPCASwLin: %g\n',out.RMSE(end));
img=showImgMask(out.alpha,out.opt.mask);
imwrite(img/max(img(:)),'FPCASwLin.png','png');
clear out004;

load(filename,'out021');
fprintf('for non skiped Ie\n');
t=1; noSkiped(:,t)=1:2000;
for i=1:3
    out=out021{1,i};
    t=t+1; noSkiped(1:length(out.cost),t)=out.cost;
    t=t+1; noSkiped(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; noSkiped(1:length(out.time),t)=out.time;
end
mincost=reshape(noSkiped(:,[2,5,8]),[],1); mincost=min(mincost(mincost>0));
noSkiped(:,[2,5,8])=noSkiped(:,[2,5,8])-mincost;
save('costRmseTime.data','noSkiped','-ascii');

t=1; skipped(:,t)=1:2000;
for i=1:3
    out=out021{2,i};
    t=t+1; skipped(1:length(out.cost),t)=out.cost;
    t=t+1; skipped(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; skipped(1:length(out.time),t)=out.time;
end
mincost=reshape(skipped(:,[2,5,8]),[],1); mincost=min(mincost(mincost>0));
skipped(:,[2,5,8])=skipped(:,[2,5,8])-mincost;
save('costRmseTime_fixedI.data','skipped','-ascii');

!find ./ \( -iname "*.png" -o -iname "*.data" \) -a -mmin -30 -exec mv -v {} /home/renliang/research/NDEReports/NDE2014_1/poster/ \;

%% runlist = 30
if(0)
    figure;
    out=out30; p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.RMSE,style1{k});
        p=p+out{k}.p;
        hold on;
    end
end

%% runlist = 31
if(0)
    figure;
    out=out31; p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.RMSE,style1{1});
        p=p+out{k}.p;
        hold on;
    end
    p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.l1Pen,style1{2});
        p=p+out{k}.p;
        hold on;
    end
    p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.cost,style1{3});
        p=p+out{k}.p;
        hold on;
    end
    p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.nonneg,style1{4});
        p=p+out{k}.p;
        hold on;
    end
end

