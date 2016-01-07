function slPoiLinEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Skyline Poisson Measurements With Linear Model Example
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
        clear -regexp '(?i)opt'
        filename = [mfilename '.mat'];
        OPT.maxItr=2e4; OPT.thresh=1e-6; OPT.debugLevel=1;
        OPT.noiseType='poisson'; OPT.matrixType='nonneg'; OPT.snr=1e7;
        OPT.contShrnk=0.1;
        m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
        a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
        aa =10.^(3:-1:-10);

        for k=1:5
            for i=1:length(m)
                OPT.m=m(i);
                [y,Phi,Phit,Psi,Psit,OPT,~,invEAAt]=loadLinear(OPT,k*100+i);
                initSig=randn(size(initSig));
                initSig=Phit(invEAAt*y)*0+1;
                fprintf('min=%d, max=%d, sum(y)=%d\n',min(y), max(y),sum(y));
                u_max=1;

                for j=5:-1:1;
                    fprintf('%s, i=%d, j=%d, k=%d\n','Skyline Poisson Linear Example',i,j,k);
                    OPT.u = a(i)*u_max*10^(j-3);

                    % % This method by Dupe etc 2009 seems not working at all
                    % dupe   {i,j,k}=Wrapper.gaussStabProxite(Phi,Phit,Psi,Psit,y,initSig,opt);
                if(i==2 && j==3 && k==1) 
                    opt=OPT; opt.thresh=1e-12; opt.maxItr=1e5;
                    spiral12{i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
                    pnpg12  {i,j,k}=Wrapper.PNPG  (Phi,Phit,Psi,Psit,y,initSig,opt);
                end

                    opt=OPT; opt.proximal='wvltFADMM';
                    fpnpg      {i,j,k}=Wrapper.PNPG     (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt=OPT; opt.alg='N83';
                    tfocs_n83_m6 {i,j,k}=Wrapper.tfocs    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt=OPT; opt.restartEvery=200;
                    tfocs_200_m6 {i,j,k}=Wrapper.tfocs    (Phi,Phit,Psi,Psit,y,initSig,opt);

                    save(filename);
                    continue
                   
                    opt=OPT;
                    pnpg      {i,j,k}=Wrapper.PNPG     (Phi,Phit,Psi,Psit,y,initSig,opt);
                    pnpgc     {i,j,k}=Wrapper.PNPGc    (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npg       {i,j,k}=Wrapper.NPG      (Phi,Phit,Psi,Psit,y,initSig,opt);
                    spiral    {i,j,k}=Wrapper.SPIRAL   (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt=OPT; opt.adaptiveStep=false;
                    pnpg_noAdp{i,j,k}=Wrapper.PNPG     (Phi,Phit,Psi,Psit,y,initSig,opt);
                    opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
                    pnpg_cumu0{i,j,k}=Wrapper.PNPG     (Phi,Phit,Psi,Psit,y,initSig,opt);
                    save(filename);
                    continue;

                end
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
            end
        end
        npgContRMSE=mean(npgContRMSE,3);

        figure;
        for i=1:length(m)
            loglog(aa,npgContRMSE(:,i),'r'); hold on;
            aaa(i)=min(npgContRMSE(:,i));
        end
        figure; semilogy(aaa); hold on; semilogy(findBest(pnpg),'ro-');

        K=5;
        fprintf('Poisson example\n');

        figure;
        semilogy(m,findBest(      pnpg,'RMSE'),'r-*' ); hold on;
        semilogy(m,findBest(       npg,'RMSE'),'c-p' );
        semilogy(m,findBest(    spiral,'RMSE'),'k-^' );
        semilogy(m,findBest(pnpg_noAdp,'RMSE'),'k*-.');
        semilogy(m,findBest(pnpg_cumu0,'RMSE'),'bs-.');
        semilogy(m,findBest(   spiral8,'RMSE'),'go-.');
        legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0','spiral8');

        figure;
        semilogy(m,findBest(      pnpg,'time'),'r-*' ); hold on;
        semilogy(m,findBest(       npg,'time'),'c-p' );
        semilogy(m,findBest(    spiral,'time'),'k-^' );
        semilogy(m,findBest(pnpg_noAdp,'time'),'k*-.');
        semilogy(m,findBest(pnpg_cumu0,'time'),'bs-.');
        semilogy(m,findBest(   spiral8,'time'),'go-.');
        legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0','spiral8');

        forSave=[findBest(      pnpg,'RMSE'),...
            findBest(       npg,'RMSE'),...
            findBest(    spiral,'RMSE'),...
            findBest(pnpg_noAdp,'RMSE'),...
            findBest(pnpg_cumu0,'RMSE'),...
            findBest(   spiral8,'RMSE'),...
            findBest(      pnpg,'time'),...
            findBest(       npg,'time'),...
            findBest(    spiral,'time'),...
            findBest(pnpg_noAdp,'time'),...
            findBest(pnpg_cumu0,'time'),...
            findBest(   spiral8,'time'),...
            m(:) ];

        save('varyMeasurementPoisson.data','forSave','-ascii');

        mIdx  = 2;
        signal=npg{1}.opt.trueAlpha;
        signal=[signal,      pnpg{mIdx,3,1}.alpha];
        signal=[signal,       npg{mIdx,3,1}.alpha];
        signal=[signal,    spiral{mIdx,4,1}.alpha];
        signal=[signal,pnpg_noAdp{mIdx,3,1}.alpha];
        signal=[signal,pnpg_cumu0{mIdx,3,1}.alpha];
        save('skylinePoisson.data','signal','-ascii');
        fprintf('      PNPG: %g%%\n',      pnpg{mIdx,3,1}.RMSE(end)*100);
        fprintf('       NPG: %g%%\n',       npg{mIdx,3,1}.RMSE(end)*100);
        fprintf('    SPIRAL: %g%%\n',    spiral{mIdx,4,1}.RMSE(end)*100);
        fprintf('PNPG_noAdp: %g%%\n',pnpg_noAdp{mIdx,3,1}.RMSE(end)*100);
        fprintf('PNPG_cumu0: %g%%\n',pnpg_cumu0{mIdx,3,1}.RMSE(end)*100);

        forSave=[]; mIdx=2; aIdx=3; fields={'cost','RMSE','time'};
        forSave=addTrace(      pnpg{mIdx,aIdx},forSave,fields);
        forSave=addTrace(       npg{mIdx,aIdx},forSave,fields);
        forSave=addTrace(    spiral{mIdx,aIdx},forSave,fields);
        forSave=addTrace(pnpg_noAdp{mIdx,aIdx},forSave,fields);
        forSave=addTrace(pnpg_cumu0{mIdx,aIdx},forSave,fields);
        forSave=addTrace(   spiral8{mIdx,aIdx},forSave,fields);
        save('cost_itr.data','forSave','-ascii');

        mincost=reshape(forSave(:,[1,4,7,10,13,16]),[],1); 
        mincost=min(mincost(mincost~=0));
        idx=(forSave(:, 1)~=0); forSave(idx, 1)=(forSave(idx, 1)-mincost);
        idx=(forSave(:, 4)~=0); forSave(idx, 4)=(forSave(idx, 4)-mincost);
        idx=(forSave(:, 7)~=0); forSave(idx, 7)=(forSave(idx, 7)-mincost);
        idx=(forSave(:,10)~=0); forSave(idx,10)=(forSave(idx,10)-mincost);
        idx=(forSave(:,13)~=0); forSave(idx,13)=(forSave(idx,13)-mincost);
        idx=(forSave(:,16)~=0); forSave(idx,16)=(forSave(idx,16)-mincost);

        figure; semilogy(forSave(:,1),'r'); hold on;
        semilogy(forSave(:,4),'g'); semilogy(forSave(:,7),'b');
        semilogy(forSave(:,10),'k'); semilogy(forSave(:,13),'c'); semilogy(forSave(:,16),'b:');
        legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0','spiral8');
        figure; semilogy(forSave(:,3),forSave(:,1),'r'); hold on;
        semilogy(forSave(:,6),forSave(:,4),'g'); semilogy(forSave(:,9),forSave(:,7),'b');
        semilogy(forSave(:,12),forSave(:,10),'k'); semilogy(forSave(:,15),forSave(:,13),'c');
        semilogy(forSave(:,18),forSave(:,16),'b:');
        legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0','spiral8');
        figure; semilogy(forSave(:,3),forSave(:,2),'r'); hold on;
        semilogy(forSave(:,6),forSave(:,5),'g'); semilogy(forSave(:,9),forSave(:,8),'b');
        semilogy(forSave(:,12),forSave(:,11),'k'); semilogy(forSave(:,15),forSave(:,14),'c');
        semilogy(forSave(:,18),forSave(:,17),'b:');
        legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0','spiral8');

        printf('\nstart to vary a to see the performance\n');

        % 3 groups are enough
        idx=[1 2 6];
        m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
        a=[1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3];
        as=1:5;
        forSave=[]; forTime=[];
        for mIdx=idx
            figure(900);
            semilogy(log10(a(mIdx))+as-3,meanOverK(      pnpg(mIdx,as,:),'RMSE'),'r-*'); hold on;
            semilogy(log10(a(mIdx))+as-3,meanOverK(       npg(mIdx,as,:),'RMSE'),'r.-');
            semilogy(log10(a(mIdx))+as-3,meanOverK(    spiral(mIdx,as,:),'RMSE'),'r-^');
            semilogy(log10(a(mIdx))+as-3,meanOverK(pnpg_noAdp(mIdx,as,:),'RMSE'),'k-s');
            semilogy(log10(a(mIdx))+as-3,meanOverK(pnpg_cumu0(mIdx,as,:),'RMSE'),'b-p');
            legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0');

            forSave=[forSave log10(a(mIdx))+as(:)-3];
            forSave=[forSave reshape(meanOverK(      pnpg(mIdx,as,:),'RMSE'),[],1)];
            forSave=[forSave reshape(meanOverK(       npg(mIdx,as,:),'RMSE'),[],1)];
            forSave=[forSave reshape(meanOverK(    spiral(mIdx,as,:),'RMSE'),[],1)];
            forSave=[forSave reshape(meanOverK(pnpg_noAdp(mIdx,as,:),'RMSE'),[],1)];
            forSave=[forSave reshape(meanOverK(pnpg_cumu0(mIdx,as,:),'RMSE'),[],1)];

            figure(901);
            semilogy(log10(a(mIdx))+as-3,meanOverK(      pnpg(mIdx,as,:),'time'),'r-*'); hold on;
            semilogy(log10(a(mIdx))+as-3,meanOverK(       npg(mIdx,as,:),'time'),'r.-');
            semilogy(log10(a(mIdx))+as-3,meanOverK(    spiral(mIdx,as,:),'time'),'r-^');
            semilogy(log10(a(mIdx))+as-3,meanOverK(pnpg_noAdp(mIdx,as,:),'time'),'g-o');
            semilogy(log10(a(mIdx))+as-3,meanOverK(pnpg_cumu0(mIdx,as,:),'time'),'k-s');
            legend('pnpg','npg','spiral','pnpg\_noAdp','pnpg\_cumu0');
            title(sprintf('mIdx=%d',mIdx));

            forTime=[forTime log10(a(mIdx))+as(:)-3];
            forTime=[forTime reshape(meanOverK(      pnpg(mIdx,as,:),'time'),[],1)];
            forTime=[forTime reshape(meanOverK(       npg(mIdx,as,:),'time'),[],1)];
            forTime=[forTime reshape(meanOverK(    spiral(mIdx,as,:),'time'),[],1)];
            forTime=[forTime reshape(meanOverK(pnpg_noAdp(mIdx,as,:),'time'),[],1)];
            forTime=[forTime reshape(meanOverK(pnpg_cumu0(mIdx,as,:),'time'),[],1)];
        end
        save('rmseVsAPoisson.data','forSave','-ascii');
        %save('timeVsAPoisson.data','forTime','-ascii');  % useless
end

end

function [a,b,c]=meanOverK(method,field)
    if(nargin==2)
        a=mean(Cell.getField(method,field),3);
    else
        a=mean(Cell.getField(method,'time'),3);
        b=mean(Cell.getField(method,'cost'),3);
        c=mean(Cell.getField(method,'RMSE'),3);
    end
end
function idx = findBestJ(method)
    rmse=meanOverK(method,'RMSE');
    [r,c]=find(rmse==repmat(min(rmse,[],2),1,5));
    [r,idx]=sort(r);
    c=c(idx);
    [r,ia]=unique(r);
    idx=c(ia);
end
function ret = findBest(method,field,idx)
    if ~exist('field','var')
        field='RMSE';
    end
    if ~exist('idx','var')
        idx=findBestJ(method);
    end
    idx=idx(:);
    data=meanOverK(method,field);
    ret=ones(size(idx))*nan;
    for i=1:min(size(data,1),length(idx))
        ret(i)=data(i,idx(i));
    end
end
function forSave=addTrace(method,forSave,fields)
    if(~exist('fields','var'))
        fields={'time','cost','RMSE'};
    end
    n=length(fields);
    for i=1:n
        data(:,i)=reshape(getfield(method,fields{i}),[],1);
    end
    forSave=appendColumns(data,forSave);
end
function forSave = appendColumns(col,forSave)
    [r,c]=size(forSave);
    forSave(1:size(col,1),c+1:c+size(col,2))=col;
end



