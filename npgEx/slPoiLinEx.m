function slPoiLinEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%
%           Skyline Poisson Linear example, no background noise
%           Vary the number of measurements, with continuation

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt'); filename = [mfilename '.mat'];
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
        opt.noiseType='poisson'; opt.matrixType='nonneg'; opt.snr=1e7;
        opt.contShrnk=0.1;
        m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
        a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
        aa =10.^(3:-1:-10);

        for k=1:5
            for i=1:length(m)
                opt.m=m(i);
                [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
                initSig=randn(size(initSig));
                initSig=Phit(invEAAt*y)*0+1;
                fprintf('min=%d, max=%d, sum(y)=%d\n',min(y), max(y),sum(y));
                u_max=1;

                for j=1:5;
                    opt.u = a(i)*u_max*10^(j-3);

                  % % This method by Dupe etc 2009 seems not working at all
                  % dupe   {i,j,k}=Wrapper.gaussStabProxite(Phi,Phit,Psi,Psit,y,initSig,opt);
                   
                    fprintf('%s, i=%d, j=%d, k=%d\n','Example_006',i,j,k);
                    if(i==2 && j==4 && k==1) 
                    spiral {i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
                    temp=opt; opt.thresh=1e-8;
                    spiral8{i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt=temp;
                end
                    continue
                    npgc   {i,j,k}=Wrapper.NPGc  (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npg    {i,j,k}=Wrapper.NPG   (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npgs   {i,j,k}=Wrapper.NPGs  (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npgsc  {i,j,k}=Wrapper.NPGsc (Phi,Phit,Psi,Psit,y,initSig,opt);
                end
                save(filename);
            end
        end

    case 'plot'
        filename = [mfilename '.mat']; load(filename);
        m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
        aa =10.^(3:-1:-10);

        K=5;
        for k=1:K
            for i=1:length(m)
                npgContRMSE(:,i,k)  = [ npgFull{i,k}.contRMSE(:); npgFull{i,k}.RMSE(end)];
                npgsContRMSE(:,i,k)  = [npgsFull{i,k}.contRMSE(:);npgsFull{i,k}.RMSE(end)];
            end
        end
        npgContRMSE=mean(npgContRMSE,3);
        npgsContRMSE=mean(npgsContRMSE,3);

        figure;
        for i=1:length(m)
            loglog(aa,npgContRMSE(:,i),'r'); hold on; loglog(aa,npgsContRMSE(:,i),'b');
            aaa(i)=min(npgContRMSE(:,i));
            bbb(i)=min(npgsContRMSE(:,i));
        end
        figure; semilogy(aaa); hold on; semilogy(bbb);

        keyboard

        K=5;
        fprintf('Poisson example\n');

        npgTime    = mean(Cell.getField(    npg(:,:,1:K),'time'),3);
        npgcTime   = mean(Cell.getField(   npgc(:,:,1:K),'time'),3);
        npgsTime   = mean(Cell.getField(   npgs(:,:,1:K),'time'),3);
        npgscTime  = mean(Cell.getField(  npgsc(:,:,1:K),'time'),3);
        spiralTime = mean(Cell.getField( spiral(:,:,1:K),'time'),3);
        spiral8Time= mean(Cell.getField(spiral8(:,:,1:K),'time'),3);

        npgCost    = mean(Cell.getField(    npg(:,:,1:K),'cost'),3);
        npgcCost   = mean(Cell.getField(   npgc(:,:,1:K),'cost'),3);
        npgsCost   = mean(Cell.getField(   npgs(:,:,1:K),'cost'),3);
        npgscCost  = mean(Cell.getField(  npgsc(:,:,1:K),'cost'),3);
        spiralCost = mean(Cell.getField( spiral(:,:,1:K),'cost'),3);
        spiral8Cost= mean(Cell.getField(spiral8(:,:,1:K),'cost'),3);

        npgRMSE    = mean(Cell.getField(    npg(:,:,1:K),'RMSE'),3);
        npgcRMSE   = mean(Cell.getField(   npgc(:,:,1:K),'RMSE'),3);
        npgsRMSE   = mean(Cell.getField(   npgs(:,:,1:K),'RMSE'),3);
        npgscRMSE  = mean(Cell.getField(  npgsc(:,:,1:K),'RMSE'),3);
        spiralRMSE = mean(Cell.getField( spiral(:,:,1:K),'RMSE'),3);
        spiral8RMSE= mean(Cell.getField(spiral8(:,:,1:K),'RMSE'),3);

        aIdx=4;
        figure;
        semilogy(m,    npgRMSE(:,aIdx),'r-*'); hold on;
        semilogy(m,   npgsRMSE(:,aIdx),'c-p');
        semilogy(m, spiralRMSE(:,aIdx),'k-^');
        semilogy(m,   npgcRMSE(:,aIdx),'k*-.');
        semilogy(m,  npgscRMSE(:,aIdx),'bs-.');
        semilogy(m,spiral8RMSE(:,aIdx),'go-.');
        legend('npg','npgs','spiral','npgc','npgsc','spiral8');

        figure;
        semilogy(m,    npgTime(:,aIdx),'r-*' ); hold on;
        semilogy(m,   npgsTime(:,aIdx),'c-p' );
        semilogy(m, spiralTime(:,aIdx),'k-^' );
        semilogy(m,   npgcTime(:,aIdx),'k*-.');
        semilogy(m,  npgscTime(:,aIdx),'bs-.');
        semilogy(m,spiral8Time(:,aIdx),'go-.');
        legend('npg','npgs','spiral','npgc','npgsc','spiral8');

        forSave=[npgTime(:,aIdx), npgsTime(:,aIdx), npgcTime(:,aIdx), npgscTime(:,aIdx), spiralTime(:,aIdx), spiral8Time(:,aIdx),...
            npgCost(:,aIdx), npgsCost(:,aIdx), npgcCost(:,aIdx), npgscCost(:,aIdx), spiralCost(:,aIdx), spiral8Cost(:,aIdx),...
            npgRMSE(:,aIdx), npgsRMSE(:,aIdx), npgcRMSE(:,aIdx), npgscRMSE(:,aIdx), spiralRMSE(:,aIdx), spiral8RMSE(:,aIdx),...
            m(:)];
        save('varyMeasurementPoisson.data','forSave','-ascii');

        temp = 2;
        signal=npg{1}.opt.trueAlpha;
        signal=[signal,    npg{temp,3,1}.alpha];
        signal=[signal,   npgc{temp,3,1}.alpha];
        signal=[signal, spiral{temp,4,1}.alpha];
        signal=[signal,   npgs{temp,3,1}.alpha];
        signal=[signal,  npgsc{temp,3,1}.alpha];
        save('skylinePoisson.data','signal','-ascii');
        fprintf('   NPG: %g%%\n',           npg{temp,3,1}.RMSE(end)*100);
        fprintf('  NPGc: %g%%\n',          npgc{temp,3,1}.RMSE(end)*100);
        fprintf('SPIRAL: %g%%\n',        spiral{temp,4,1}.RMSE(end)*100);
        fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{temp,3,1}.RMSE(end)*100,rmseTruncate( npgs{temp,3,1})*100);
        fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{temp,3,1}.RMSE(end)*100,rmseTruncate(npgsc{temp,3,1})*100);

        forSave=[]; t=0; mIdx=2;
        out=  npgc{mIdx,aIdx,1};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=   npg{mIdx,aIdx,1};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=spiral8{mIdx,aIdx,1};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=   npgs{mIdx,aIdx,1};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=  npgsc{mIdx,aIdx,1};
        t=t+1; forSave(1:length(out.cost),t)=out.cost;
        t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
        t=t+1; forSave(1:length(out.time),t)=out.time;
        out=spiral8{mIdx,aIdx,1}; keyboard
        q=(1:max(find(spiral8{mIdx,aIdx}.cost(:)>=spiral{mIdx,aIdx}.cost(end))))';
        t=t+1; forSave(1:length(out.cost(q)),t)=out.cost(q);
        t=t+1; forSave(1:length(out.RMSE(q)),t)=out.RMSE(q);
        t=t+1; forSave(1:length(out.time(q)),t)=out.time(q);

        save('cost_itr.data','forSave','-ascii');
        mincost=reshape(forSave(:,[1,4,7,10,13,16]),[],1); 
        mincost=min(mincost(mincost~=0));
        idx=(forSave(:, 1)~=0); forSave(idx, 1)=(forSave(idx, 1)-mincost);
        idx=(forSave(:, 4)~=0); forSave(idx, 4)=(forSave(idx, 4)-mincost);
        idx=(forSave(:, 7)~=0); forSave(idx, 7)=(forSave(idx, 7)-mincost);
        idx=(forSave(:,10)~=0); forSave(idx,10)=(forSave(idx,10)-mincost);
        idx=(forSave(:,13)~=0); forSave(idx,13)=(forSave(idx,13)-mincost);

        figure; semilogy(forSave(:,1),'r'); hold on;
        semilogy(forSave(:,4),'g'); semilogy(forSave(:,7),'b');
        semilogy(forSave(:,10),'k'); semilogy(forSave(:,13),'c');
        legend('npgc','npg','spiral8','npgs','npgsc');
        figure; loglog(forSave(:,3),forSave(:,1),'r'); hold on;
        semilogy(forSave(:,6),forSave(:,4),'g'); semilogy(forSave(:,9),forSave(:,7),'b');
        semilogy(forSave(:,12),forSave(:,10),'k'); semilogy(forSave(:,15),forSave(:,13),'c');
        legend('npgc','npg','spiral8','npgs','npgsc');
        figure; loglog(forSave(:,3),forSave(:,2),'r'); hold on;
        semilogy(forSave(:,6),forSave(:,5),'g'); semilogy(forSave(:,9),forSave(:,8),'b');
        semilogy(forSave(:,12),forSave(:,11),'k'); semilogy(forSave(:,15),forSave(:,14),'c');
        legend('npgc','npg','spiral8','npgs','npgsc');

        keyboard

        printf('\nstart to vary a to see the performance\n');

        % 3 groups are enough
        idx=[1 2 6];
        m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
        a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
        as=1:5;
        forSave=[]; forTime=[];
        for mIdx=idx
            figure(900);
            semilogy(log10(a(mIdx))+as-3,    npgRMSE(mIdx,as),'r-*'); hold on;
            semilogy(log10(a(mIdx))+as-3,   npgcRMSE(mIdx,as),'r.-');
            semilogy(log10(a(mIdx))+as-3, spiralRMSE(mIdx,as),'r-^');
            semilogy(log10(a(mIdx))+as-3,spiral8RMSE(mIdx,as),'g-o');
            semilogy(log10(a(mIdx))+as-3,   npgsRMSE(mIdx,as),'k-s');
            semilogy(log10(a(mIdx))+as-3,  npgscRMSE(mIdx,as),'b-p');
            legend('npg','npgc','spiral','spiral8','npgs','npgsc');

            forSave=[forSave log10(a(mIdx))+as(:)-3];
            forSave=[forSave reshape(    npgRMSE(mIdx,as),[],1)];
            forSave=[forSave reshape(   npgcRMSE(mIdx,as),[],1)];
            forSave=[forSave reshape( spiralRMSE(mIdx,as),[],1)];
            forSave=[forSave reshape(spiral8RMSE(mIdx,as),[],1)];
            forSave=[forSave reshape(   npgsRMSE(mIdx,as),[],1)];
            forSave=[forSave reshape(  npgscRMSE(mIdx,as),[],1)];

            figure(901);
            semilogy(log10(a(mIdx))+as-3,    npgTime(mIdx,as),'r-*'); hold on;
            semilogy(log10(a(mIdx))+as-3,   npgcTime(mIdx,as),'r.-');
            semilogy(log10(a(mIdx))+as-3, spiralTime(mIdx,as),'r-^');
            semilogy(log10(a(mIdx))+as-3,spiral8Time(mIdx,as),'g-o');
            semilogy(log10(a(mIdx))+as-3,   npgsTime(mIdx,as),'k-s');
            semilogy(log10(a(mIdx))+as-3,  npgscTime(mIdx,as),'b-p');
            legend('npg','npgc','spiral','spiral8','npgs','npgsc');
            title(sprintf('mIdx=%d',mIdx));

            forTime=[forTime log10(a(mIdx))+as(:)-3];
            forTime=[forTime reshape(    npgTime(mIdx,as),[],1)];
            forTime=[forTime reshape(   npgcTime(mIdx,as),[],1)];
            forTime=[forTime reshape( spiralTime(mIdx,as),[],1)];
            forTime=[forTime reshape(spiral8Time(mIdx,as),[],1)];
            forTime=[forTime reshape(   npgsTime(mIdx,as),[],1)];
            forTime=[forTime reshape(  npgscTime(mIdx,as),[],1)];
        end
        save('rmseVsAPoisson.data','forSave','-ascii');
        save('timeVsAPoisson.data','forTime','-ascii');

        system(['mv rmseVsAPoisson.data timeVsAPoisson.data varyMeasurementPoisson.data cost_itr.data ' paperDir]);
        system(['mv skylinePoisson.data ' paperDir]);

    case 'test'  % run to get good choice of a
        filename = [mfilename '_test.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt');
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
        opt.noiseType='poisson'; opt.matrixType='nonneg'; opt.snr=1e7;
        m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
        a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
        aa =10.^(3:-1:-10);

        for k=1:5
            for i=1:length(m)
                opt.m=m(i);
                [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
                initSig=randn(size(initSig));
                initSig=Phit(invEAAt*y)*0+1;
                fprintf('min=%d, max=%d, sum(y)=%d\n',min(y), max(y),sum(y));
                u_max=1;

                opt.fullcont=true;
                opt.u=aa*u_max;
                npgFull {i,k}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgFull{i,k};
                fprintf('i=%d, good a = 1e%g\n',i,max(log10(aa(out.contRMSE==min(out.contRMSE)))));
                npgsFull{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt); out=npgsFull{i,k};
                fprintf('i=%d, good a = 1e%g\n',i,max(log10(aa(out.contRMSE==min(out.contRMSE)))));
                opt.fullcont=false;
                save(filename);
            end
        end
end

end
