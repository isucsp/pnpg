function runAsilomar2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration

if(nargin==0) runList = [0]; elseif(isempty(runList)) return; end
paperDir = '~/research/myPaper/asilomar2014/';

if(any(runList==002))
    filename = [mfilename '_002.mat']; load(filename);

    m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    u = [1e-3,1e-3,1e-4,1e-4,1e-5,1e-5,1e-6,1e-6,1e-6];
    idx=2:2:7;
    K = 1;

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    npgscTime = mean(Cell.getField( npgsc(:,:,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,:,1:K),'time'),3);
    fpcasTime = mean(Cell.getField( fpcas(:,:,1:K),'cpu' ),3);
    fpcTime   = mean(Cell.getField(   fpc(:,:,1:K),'time' ),3);
    fistaTime = mean(Cell.getField( fista(:,:,1:K),'time'),3);
    sparsaTime= mean(Cell.getField(sparsa(:,:,1:K),'time'),3);
    sparsnTime= mean(Cell.getField(sparsn(:,:,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    npgscCost = mean(Cell.getField( npgsc(:,:,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,:,1:K),'cost'),3);
    fpcasCost = mean(Cell.getField( fpcas(:,:,1:K),'f'   ),3);
    fpcCost   = mean(Cell.getField(   fpc(:,:,1:K),'cost' ),3);
    fistaCost = mean(Cell.getField( fista(:,:,1:K),'cost'),3);
    sparsaCost= mean(Cell.getField(sparsa(:,:,1:K),'cost'),3);
    sparsnCost= mean(Cell.getField(sparsn(:,:,1:K),'cost'),3);

    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,:,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,:,1:K),'reconerror'),3);
    fpcasRMSE = mean(Cell.getField( fpcas(:,:,1:K),'RMSE'),3);
    fpcRMSE   = mean(Cell.getField(   fpc(:,:,1:K),'RMSE' ),3);
    fistaRMSE = mean(Cell.getField( fista(:,:,1:K),'RMSE'),3);
    sparsaRMSE= mean(Cell.getField(sparsa(:,:,1:K),'RMSE'),3);
    sparsnRMSE= mean(Cell.getField(sparsn(:,:,1:K),'RMSE'),3);

    npgnasTime = mean(Cell.getField(   npg_nads(:,:,1:K),'time'),3);
    npgcnasTime= mean(Cell.getField(  npgc_nads(:,:,1:K),'time'),3);
    npgnasCost = mean(Cell.getField(   npg_nads(:,:,1:K),'cost'),3);
    npgcnasCost= mean(Cell.getField(  npgc_nads(:,:,1:K),'cost'),3);
    npgnasRMSE = mean(Cell.getField(   npg_nads(:,:,1:K),'RMSE'),3);
    npgcnasRMSE= mean(Cell.getField(  npgc_nads(:,:,1:K),'RMSE'),3);

    for i=1:length(m)
        temp=[];
        for k=1:K
            temp(:,k)=gnet{i,k}.RMSE(:);
        end
        %gnetRMSE(i,:)=mean(temp,2)';
    end

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c4]=find(spiralRMSE== repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c8]=find(sparsnRMSE== repmat(min(sparsnRMSE,[],2),1,5)); [r,idx8]=sort(r);

    [r,c3]=find(  npgsRMSE== repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r);
    [r,c7]=find(sparsaRMSE== repmat(min(sparsaRMSE,[],2),1,5)); [r,idx7]=sort(r);
    disp([c1(idx1), c2(idx2), c4(idx4), c8(idx8) zeros(9,1) c3(idx3), c5(idx5), c6(idx6), c7(idx7) ]);
    keyboard
    uNonneg=[3 3 3 3 4 4 4 4 3];
       uNeg=[4 4 4 4 4 4 4 4 3];
    figure;
    semilogy(m,   npgRMSE((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcRMSE((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsRMSE((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralRMSE((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m, fpcasRMSE((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaRMSE((c6(idx6)-1)*9+(1:9)'),'b-.');
    semilogy(m,sparsaRMSE((c7(idx7)-1)*9+(1:9)'),'y-p');
    semilogy(m,sparsnRMSE((c8(idx8)-1)*9+(1:9)'),'r-x');
    legend('npg','npgc','npgs','spiral','fpcas','fista','sparas','sparsa');
    figure;
    semilogy(m,   npgTime((c1(idx1)-1)*9+(1:9)'),'r-*'); hold on;
    semilogy(m,  npgcTime((c2(idx2)-1)*9+(1:9)'),'c-p');
    semilogy(m,  npgsTime((c3(idx3)-1)*9+(1:9)'),'k-s');
    semilogy(m,spiralTime((c4(idx4)-1)*9+(1:9)'),'k-^');
    semilogy(m, fpcasTime((c5(idx5)-1)*9+(1:9)'),'g-o');
    semilogy(m, fistaTime((c6(idx6)-1)*9+(1:9)'),'b-.');
    semilogy(m,sparsaTime((c7(idx7)-1)*9+(1:9)'),'y-p');
    semilogy(m,sparsnTime((c8(idx8)-1)*9+(1:9)'),'r-x');
    legend('npg','npgc','npgs','spiral','fpcas','fista','sparas','sparsa');
    figure;
    semilogy(m,   npgRMSE((uNonneg-1)*9+(1:9)),'r-*'); hold on;
    semilogy(m,  npgcRMSE((uNonneg-1)*9+(1:9)),'c-p');
    semilogy(m,spiralRMSE((uNonneg-1)*9+(1:9)),'k-^');
    semilogy(m,sparsnRMSE((uNonneg-1)*9+(1:9)),'r-x');
    semilogy(m,  npgsRMSE((uNeg   -1)*9+(1:9)),'k-s');
    semilogy(m, fpcasRMSE((uNeg   -1)*9+(1:9)),'g-o');
    semilogy(m, fistaRMSE((uNeg   -1)*9+(1:9)),'b-.');
    semilogy(m,sparsaRMSE((uNeg   -1)*9+(1:9)),'y-p');
    %semilogy(m,  gnetRMSE((uNeg   -1)*9+(1:9)),'r:>');
    legend('npg','npgc','spiral','sparsa','npgs','fpcas','fista','sparas');
    figure;
    semilogy(m,   npgTime((uNonneg-1)*9+(1:9)),'r-*'); hold on;
    semilogy(m,  npgcTime((uNonneg-1)*9+(1:9)),'c-p');
    semilogy(m,spiralTime((uNonneg-1)*9+(1:9)),'k-^');
    semilogy(m,sparsnTime((uNonneg-1)*9+(1:9)),'r-x');
    semilogy(m,  npgsTime((uNeg   -1)*9+(1:9)),'k-s');
    semilogy(m, fpcasTime((uNeg   -1)*9+(1:9)),'g-o');
    semilogy(m, fistaTime((uNeg   -1)*9+(1:9)),'b-.');
    semilogy(m,sparsaTime((uNeg   -1)*9+(1:9)),'y-p');
    legend('npg','npgc','spiral','sparsa','npgs','fpcas','fista','sparas');

    f=fopen('selectedTime.data','w');
    for mIdx=1:length(m)
        fprintf(f,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%s\t%s\n',...
                npgTime(mIdx,uNonneg(mIdx)), ...
               npgcTime(mIdx,uNonneg(mIdx)), ...
             spiralTime(mIdx,uNonneg(mIdx)), ...
             sparsnTime(mIdx,uNonneg(mIdx)), ...
               npgsTime(mIdx,uNeg   (mIdx)), ...
              fpcasTime(mIdx,uNeg   (mIdx)), ...
              fistaTime(mIdx,uNeg   (mIdx)), ...
             sparsaTime(mIdx,uNeg   (mIdx)), ...
            m(mIdx),num2str(log10(u((mIdx)))+uNonneg(mIdx)-3), num2str(log10(u(mIdx))+uNeg(mIdx)-3));
    end
    fclose(f);

    as=1:5;
    forSave=[]; forTime=[];
    for mIdx=idx
        figure(900);
        semilogy(log10(u(mIdx))+as-3,    npgRMSE(mIdx,as),'r-*'); hold on;
        semilogy(log10(u(mIdx))+as-3,   npgcRMSE(mIdx,as),'r.-');
        semilogy(log10(u(mIdx))+as-3, sparsnRMSE(mIdx,as),'r-s');
        semilogy(log10(u(mIdx))+as-3, spiralRMSE(mIdx,as),'r-^');
        semilogy(log10(u(mIdx))+as-3,   npgsRMSE(mIdx,as),'k-s');
        semilogy(log10(u(mIdx))+as-3,  fpcasRMSE(mIdx,as),'g-o');
        semilogy(log10(u(mIdx))+as-3,  fistaRMSE(mIdx,as),'g-.');
        semilogy(log10(u(mIdx))+as-3, sparsaRMSE(mIdx,as),'g->');
        semilogy(log10(u(mIdx))+as-3,  npgscRMSE(mIdx,as),'g-*');
        semilogy(log10(u(mIdx))+as-3,    fpcRMSE(mIdx,as),'g:p');
        semilogy(gnet{mIdx,1}.a,   gnet{mIdx,1}.RMSE(:),'r:>');

        forSave=[forSave log10(u(mIdx))+as(:)-3];
        forSave=[forSave reshape(   npgRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(  npgcRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(sparsnRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(spiralRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(  npgsRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape( fpcasRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape( fistaRMSE(mIdx,as),[],1)];
        forSave=[forSave reshape(sparsaRMSE(mIdx,as),[],1)];

        figure;
        semilogy(log10(u(mIdx))+as-3,    npgTime(mIdx,as),'r-*'); hold on;
        semilogy(log10(u(mIdx))+as-3,   npgcTime(mIdx,as),'r.-');
        semilogy(log10(u(mIdx))+as-3, sparsnTime(mIdx,as),'r-s');
        semilogy(log10(u(mIdx))+as-3, spiralTime(mIdx,as),'r-^');
        semilogy(log10(u(mIdx))+as-3,   npgsTime(mIdx,as),'k-s');
        semilogy(log10(u(mIdx))+as-3,  fpcasTime(mIdx,as),'g-o');
        semilogy(log10(u(mIdx))+as-3,  fistaTime(mIdx,as),'g-.');
        semilogy(log10(u(mIdx))+as-3, sparsaTime(mIdx,as),'g->');
        semilogy(log10(u(mIdx))+as-3,  npgscTime(mIdx,as),'g-*');
        semilogy(log10(u(mIdx))+as-3,    fpcTime(mIdx,as),'g:p');
        legend('npg','npgc','sparsn','spiral','npgs','fpcas','fista','sparas');
        title(sprintf('mIdx=%d',mIdx));

        forTime=[forTime log10(u(mIdx))+as(:)-3];
        forTime=[forTime reshape(   npgTime(mIdx,as),[],1)];
        forTime=[forTime reshape(  npgcTime(mIdx,as),[],1)];
        forTime=[forTime reshape(sparsnTime(mIdx,as),[],1)];
        forTime=[forTime reshape(spiralTime(mIdx,as),[],1)];
        forTime=[forTime reshape(  npgsTime(mIdx,as),[],1)];
        forTime=[forTime reshape( fpcasTime(mIdx,as),[],1)];
        forTime=[forTime reshape( fistaTime(mIdx,as),[],1)];
        forTime=[forTime reshape(sparsaTime(mIdx,as),[],1)];
    end
    figure(900); 
    legend('npg','npgc','sparsn','spiral','npgs','fpcas','fista','sparas','npgsc','fpc','glmnet');
    save('rmseVsA.data','forSave','-ascii');
    save('timeVsA.data','forTime','-ascii');

    keyboard

    mIdx=6; as=gEle(c2(idx2),mIdx); forSave=[]; t=0; q=(1:max(find(sparsn12{mIdx,as}.cost(:)>=sparsn{mIdx,as}.cost(end))))';
    t=t+1; temp=      npg{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      npg{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.RMSE(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.time(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.cost(q);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= sparsn12{mIdx,as}.difAlpha(q); forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=     npgs{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=    npgsc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=   sparsa{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   sparsa{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   sparsa{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   sparsa{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=      fpc{mIdx,as}.difAlpha(:); forSave(1:length(temp),t)=temp;
    save('traceLinGauss.data','forSave','-ascii');

    mc=forSave(:,[15,19,23,27]); mc = min(mc(mc(:)>0));
    figure;
    semilogy(forSave(:,14),forSave(:,15)-mc,'r-'); hold on;
    semilogy(forSave(:,18),forSave(:,19)-mc,'r-.');
    semilogy(forSave(:,22),forSave(:,23)-mc,'c--');
    semilogy(forSave(:,26),forSave(:,27)-mc,'b:'); 

    mIdx=6; as=gEle(c2(idx2),mIdx); forSave=[]; t=0;
    t=t+1; temp=  npgc{mIdx,as}.RMSE(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.time(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.cost(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.difAlpha(:);  forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.uRecord(:,2); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgc{mIdx,as}.contThresh(:);forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.RMSE(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.time(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.cost(:);      forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.difAlpha(:);  forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.uRecord(:,2); forSave(1:length(temp),t)=temp;
    t=t+1; temp= npgsc{mIdx,as}.contThresh(:);forSave(1:length(temp),t)=temp;
    save('continuation.data','forSave','-ascii');

    keyboard

    temp = 4;
    signal=npg{1}.opt.trueAlpha;
    signal=[signal,    npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,   npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,   npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  fista{gEle((c6(idx6)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, sparsa{gEle((c7(idx7)-1)*9+(1:9)',temp)}.alpha];
    save('skyline.data','signal','-ascii');
    figure; plot(signal(:,2)); hold on; plot(signal(:,1),'r'); title('NPG');
    figure; plot(signal(:,4)); hold on; plot(signal(:,1),'r'); title('NPGs');
    figure; plot(signal(:,6)); hold on; plot(signal(:,1),'r'); title('FPCas');
    fprintf('npgRec RMSE: %g%%\n',npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100);
    fprintf('npgsRec RMSE: %g%%\n',npgs{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100);
    fprintf('fistaRec RMSE: %g%%\n',fista{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100);

    M=length(m);
    str=        '$m$            ';              for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%10d',str,m(i)); end; end;
    str=sprintf('%s\\\\\\hline',str);
    str=sprintf('%s\\\\\nNPG            ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,   npg{gEle((c1(idx1)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nNPG$_\\text{C}$ ',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,  npgc{gEle((c2(idx2)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nNPG$_\\text{S}$ ',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,  npgs{gEle((c3(idx3)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nSPIRAL         ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,spiral{gEle((c4(idx4)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nFPC$_\\text{AS}$',str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str, fpcas{gEle((c5(idx5)-1)*9+(1:9)',i)}.f   (end));end; end;
    str=sprintf('%s\nFISTA          ', str);    for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str, fista{gEle((c6(idx6)-1)*9+(1:9)',i)}.cost(end));end; end;
    str=sprintf('%s\\\\\nSpaRSA         ', str);for i=1:M;if(mod(m(i),100)==0);str=sprintf('%s&%-10.4g',str,spiral{gEle((c7(idx7)-1)*9+(1:9)',i)}.cost(end));end; end;
    file=fopen('varyMeasurementTable.tex','w'); fprintf(file,'%s',str); fclose(file);

    % figure;
    % for i=1:M;
    %     semilogy(npgs{i,idx3,1}.stepSize); hold on; semilogy(fista{i,idx6,1}.stepSize,'r:');
    %     semilogy([1,length(fista{i,idx6,1}.RMSE)],ones(1,2)*1/fista{i,idx6,1}.opt.L,'k-.');
    %     hold off;
    %     pause;
    % end
      npgItr=[];   
     npgcItr=[];
     npgsItr=[];
   spiralItr=[];
    fpcasItr=[];
    fistaItr=[];
   sparsaItr=[];
   sparsnItr=[];

    for i=1:K
        temp=   npg(:,:,i); temp=temp((c1(idx1)-1)*9+(1:9)');    npgItr=[   npgItr,showResult(temp,2,'p'   )];
        temp=  npgc(:,:,i); temp=temp((c2(idx2)-1)*9+(1:9)');   npgcItr=[  npgcItr,showResult(temp,2,'p'   )];
        temp=  npgs(:,:,i); temp=temp((c3(idx3)-1)*9+(1:9)');   npgsItr=[  npgsItr,showResult(temp,2,'p'   )];
        temp=spiral(:,:,i); temp=temp((c4(idx4)-1)*9+(1:9)'); spiralItr=[spiralItr,showResult(temp,2,'p'   )];
        temp= fpcas(:,:,i); temp=temp((c5(idx5)-1)*9+(1:9)');  fpcasItr=[ fpcasItr,showResult(temp,2,'itr' )];
        temp= fista(:,:,i); temp=temp((c6(idx6)-1)*9+(1:9)');  fistaItr=[ fistaItr,showResult(temp,2,'p'   )];
        temp=sparsa(:,:,i); temp=temp((c7(idx7)-1)*9+(1:9)'); sparsaItr=[sparsaItr,showResult(temp,3,'RMSE')];
        temp=sparsn(:,:,i); temp=temp((c8(idx8)-1)*9+(1:9)'); sparsnItr=[sparsnItr,showResult(temp,3,'RMSE')];
    end

    forSave=[];
    forSave=[forSave,    npgTime((c1(idx1)-1)*9+(1:9)')];
    forSave=[forSave,   npgcTime((c2(idx2)-1)*9+(1:9)')];
    forSave=[forSave,   npgsTime((c3(idx3)-1)*9+(1:9)')];
    forSave=[forSave, spiralTime((c4(idx4)-1)*9+(1:9)')];
    forSave=[forSave,  fpcasTime((c5(idx5)-1)*9+(1:9)')];
    forSave=[forSave,  fistaTime((c6(idx6)-1)*9+(1:9)')];

    forSave=[forSave,    npgCost((c1(idx1)-1)*9+(1:9)')];
    forSave=[forSave,   npgcCost((c2(idx2)-1)*9+(1:9)')];
    forSave=[forSave,   npgsCost((c3(idx3)-1)*9+(1:9)')];
    forSave=[forSave, spiralCost((c4(idx4)-1)*9+(1:9)')];
    forSave=[forSave,  fpcasCost((c5(idx5)-1)*9+(1:9)')];
    forSave=[forSave,  fistaCost((c6(idx6)-1)*9+(1:9)')];

    forSave=[forSave,    npgRMSE((c1(idx1)-1)*9+(1:9)')];
    forSave=[forSave,   npgcRMSE((c2(idx2)-1)*9+(1:9)')];
    forSave=[forSave,   npgsRMSE((c3(idx3)-1)*9+(1:9)')];
    forSave=[forSave, spiralRMSE((c4(idx4)-1)*9+(1:9)')];
    forSave=[forSave,  fpcasRMSE((c5(idx5)-1)*9+(1:9)')];
    forSave=[forSave,  fistaRMSE((c6(idx6)-1)*9+(1:9)')];
    forSave=[forSave, m(:)];
    forSave=[forSave, sparsaTime((c7(idx7)-1)*9+(1:9)')];
    forSave=[forSave, sparsaCost((c7(idx7)-1)*9+(1:9)')];
    forSave=[forSave, sparsaRMSE((c7(idx7)-1)*9+(1:9)')];
    forSave=[forSave, sparsnTime((c8(idx8)-1)*9+(1:9)')];
    forSave=[forSave, sparsnCost((c8(idx8)-1)*9+(1:9)')];
    forSave=[forSave, sparsnRMSE((c8(idx8)-1)*9+(1:9)')];
    save('varyMeasurement.data','forSave','-ascii');

    forSave=m(:);
    forSave=[forSave,    npgTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,   npgcTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave, spiralTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave, sparsnTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,   npgsTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,  fpcasTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave,  fistaTime((uNonneg-1)*9+(1:9))'];
    forSave=[forSave, sparsaTime((uNonneg-1)*9+(1:9))'];
    save('varyMeasurementTime.data','forSave','-ascii');

    keyboard

    mIdx=2; experi=1; forSave=[]; t=0;
    npgsT=npgsT(:,:,experi); npgsn20T=npgs(:,:,experi); fistaT=fista(:,:,experi); fistalT=fistal(:,:,experi); fistalT{9,6}=[];
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c6(idx6),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c6(idx6),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c6(idx6),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=   npgsT{mIdx,gEle(c3(idx3),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=npgsn20T{mIdx,gEle(c3(idx3),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp=  fistaT{mIdx,gEle(c6(idx6),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c6(idx6),mIdx)}.cost(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c6(idx6),mIdx)}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c6(idx6),mIdx)}.RMSE(:);     forSave(1:length(temp),t)=temp;
    t=t+1; temp= fistalT{mIdx,gEle(c6(idx6),mIdx)}.time(:);     forSave(1:length(temp),t)=temp;
    disp([ c3(idx3) c6(idx6)]);
    disp([   npgsT{mIdx,gEle(c3(idx3),mIdx)}.p;...
           fistalT{mIdx,gEle(c3(idx3),mIdx)}.p;...
            fistaT{mIdx,gEle(c6(idx6),mIdx)}.p]);
    disp([   npgsT{mIdx,gEle(c3(idx3),mIdx)}.time(end); ...
           fistalT{mIdx,gEle(c3(idx3),mIdx)}.time(end); ...
            fistaT{mIdx,gEle(c6(idx6),mIdx)}.time(end)]);
    temp=forSave(:,10:13); temp=temp(:); temp=temp(temp>0); temp=min(temp); forSave(:,10:13)=forSave(:,10:13)-temp;
    save('stepSizeLin.data','forSave','-ascii');
    figure(1); hold off; semilogy(forSave(:,7),forSave(:,4),'r'); hold on;
    semilogy(forSave(:,8),forSave(:,5),'g');
    semilogy(forSave(:,9),forSave(:,6),'b');
    figure(2); hold off; semilogy(forSave(:,7),forSave(:,10),'r'); hold on;
    semilogy(forSave(:,8),forSave(:,11),'g');
    semilogy(forSave(:,9),forSave(:,12),'b');
    semilogy(forSave(:,16),forSave(:,13),'c');
    keyboard

    system(['mv continuation.data traceLinGauss.data selectedTime.data timeVsA.data rmseVsA.data stepSizeLin.data varyMeasurement.data varyMeasurementTime.data skyline.data varyMeasurementTable.tex ' paperDir]);
    disp('done');
end

if(any(runList==005))
    filename = [mfilename '_005.mat']; load(filename);
    snr=[  10,  50, 100, 200, 500, 1e3, 1e4, 1e5, 1e6, 1e7]';
    u  =[1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3,1e-4,1e-4,1e-5];

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,:,1:K),'time'),3);
    fpcasTime = mean(Cell.getField( fpcas(:,:,1:K),'cpu' ),3);
    fistaTime = mean(Cell.getField( fista(:,:,1:K),'time'),3);
    sparsaTime= mean(Cell.getField(sparsa(:,:,1:K),'time'),3);
    sparsnTime= mean(Cell.getField(sparsn(:,:,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,:,1:K),'cost'),3);
    fpcasCost = mean(Cell.getField( fpcas(:,:,1:K),'f'   ),3);
    fistaCost = mean(Cell.getField( fista(:,:,1:K),'cost'),3);
    sparsaCost= mean(Cell.getField(sparsa(:,:,1:K),'cost'),3);
    sparsnCost= mean(Cell.getField(sparsn(:,:,1:K),'cost'),3);

    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,:,1:K),'reconerror'),3);
    fpcasRMSE = mean(Cell.getField( fpcas(:,:,1:K),'RMSE'),3);
    fistaRMSE = mean(Cell.getField( fista(:,:,1:K),'RMSE'),3);
    sparsaRMSE= mean(Cell.getField(sparsa(:,:,1:K),'RMSE'),3);
    sparsnRMSE= mean(Cell.getField(sparsn(:,:,1:K),'RMSE'),3);

    [r,c1]=find(   npgRMSE== repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgcRMSE== repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c4]=find(spiralRMSE== repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c8]=find(sparsnRMSE== repmat(min(sparsnRMSE,[],2),1,5)); [r,idx8]=sort(r);

    [r,c3]=find(  npgsRMSE== repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c5]=find( fpcasRMSE== repmat(min( fpcasRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find( fistaRMSE== repmat(min( fistaRMSE,[],2),1,5)); [r,idx6]=sort(r);
    [r,c7]=find(sparsaRMSE== repmat(min(sparsaRMSE,[],2),1,5)); [r,idx7]=sort(r);
    disp([c1(idx1), c2(idx2), c4(idx4), c8(idx8) zeros(10,1) c3(idx3), c5(idx5), c6(idx6), c7(idx7) ]);

    figure;
    loglog(snr,   npgRMSE((c1(idx1)-1)*10+(1:10)'),'r-*'); hold on;
    loglog(snr,  npgcRMSE((c2(idx2)-1)*10+(1:10)'),'c-p');
    loglog(snr,  npgsRMSE((c3(idx3)-1)*10+(1:10)'),'k-s');
    loglog(snr,spiralRMSE((c4(idx4)-1)*10+(1:10)'),'k-^');
    loglog(snr, fpcasRMSE((c5(idx5)-1)*10+(1:10)'),'g-o');
    loglog(snr, fistaRMSE((c6(idx6)-1)*10+(1:10)'),'b-.');
    loglog(snr,sparsaRMSE((c7(idx7)-1)*10+(1:10)'),'y-p');
    loglog(snr,sparsnRMSE((c8(idx8)-1)*10+(1:10)'),'r-x');
    legend('npg','npgc','npgs','spiral','fpcas','fista','sparas','sparsn');
    figure;
    loglog(snr,   npgTime((c1(idx1)-1)*10+(1:10)'),'r-*'); hold on;
    loglog(snr,  npgcTime((c2(idx2)-1)*10+(1:10)'),'c-p');
    loglog(snr,  npgsTime((c3(idx3)-1)*10+(1:10)'),'k-s');
    loglog(snr,spiralTime((c4(idx4)-1)*10+(1:10)'),'k-^');
    loglog(snr, fpcasTime((c5(idx5)-1)*10+(1:10)'),'g-o');
    loglog(snr, fistaTime((c6(idx6)-1)*10+(1:10)'),'b-.');
    loglog(snr,sparsaTime((c7(idx7)-1)*10+(1:10)'),'y-p');
    loglog(snr,sparsnTime((c8(idx8)-1)*10+(1:10)'),'r-x');
    legend('npg','npgc','npgs','spiral','fpcas','fista','sparas','sparsn');

    M=length(snr);
    str=        'SNR            ';              for i=1:M; str=sprintf('%s&%10d',str,snr(i)); end;
    str=sprintf('%s\\\\\\hline',str);
    str=sprintf('%s\\\\\nNPG            ', str);for i=1:M;str=sprintf('%s&%-10.4g',str,   npg{gEle((c1(idx1)-1)*10+(1:10)',i)}.cost(end));end;
    str=sprintf('%s\\\\\nNPG$_\\text{C}$ ',str);for i=1:M;str=sprintf('%s&%-10.4g',str,  npgc{gEle((c2(idx2)-1)*10+(1:10)',i)}.cost(end));end;
    str=sprintf('%s\\\\\nNPG$_\\text{S}$ ',str);for i=1:M;str=sprintf('%s&%-10.4g',str,  npgs{gEle((c3(idx3)-1)*10+(1:10)',i)}.cost(end));end;
    str=sprintf('%s\\\\\nSPIRAL         ', str);for i=1:M;str=sprintf('%s&%-10.4g',str,spiral{gEle((c4(idx4)-1)*10+(1:10)',i)}.cost(end));end;
    str=sprintf('%s\\\\\nFPC$_\\text{AS}$',str);for i=1:M;str=sprintf('%s&%-10.4g',str, fpcas{gEle((c5(idx5)-1)*10+(1:10)',i)}.f   (end));end;
    str=sprintf('%s\nFISTA          ', str);    for i=1:M;str=sprintf('%s&%-10.4g',str, fista{gEle((c6(idx6)-1)*10+(1:10)',i)}.cost(end));end;
    str=sprintf('%s\\\\\nSpaRSA         ', str);for i=1:M;str=sprintf('%s&%-10.4g',str,spiral{gEle((c7(idx7)-1)*10+(1:10)',i)}.cost(end));end;
    file=fopen('varySNRTable.tex','w'); fprintf(file,'%s',str); fclose(file);

      npgItr=[];   
     npgcItr=[];
     npgsItr=[];
   spiralItr=[];
    fpcasItr=[];
    fistaItr=[];
   sparsaItr=[];

    for i=1:K
        temp=   npg(:,:,i); temp=temp((c1(idx1)-1)*10+(1:10)');    npgItr=[   npgItr,showResult(temp,2,'p'   )];
        temp=  npgc(:,:,i); temp=temp((c2(idx2)-1)*10+(1:10)');   npgcItr=[  npgcItr,showResult(temp,2,'p'   )];
        temp=  npgs(:,:,i); temp=temp((c3(idx3)-1)*10+(1:10)');   npgsItr=[  npgsItr,showResult(temp,2,'p'   )];
        temp=spiral(:,:,i); temp=temp((c4(idx4)-1)*10+(1:10)'); spiralItr=[spiralItr,showResult(temp,2,'p'   )];
        temp= fpcas(:,:,i); temp=temp((c5(idx5)-1)*10+(1:10)');  fpcasItr=[ fpcasItr,showResult(temp,2,'itr' )];
        temp= fista(:,:,i); temp=temp((c6(idx6)-1)*10+(1:10)');  fistaItr=[ fistaItr,showResult(temp,2,'p'   )];
        temp=sparsa(:,:,i); temp=temp((c7(idx7)-1)*10+(1:10)'); sparsaItr=[sparsaItr,showResult(temp,3,'RMSE')];
    end
    keyboard

    forSave=[];
    forSave=[forSave,    npgTime((c1(idx1)-1)*10+(1:10)')];
    forSave=[forSave,   npgcTime((c2(idx2)-1)*10+(1:10)')];
    forSave=[forSave,   npgsTime((c3(idx3)-1)*10+(1:10)')];
    forSave=[forSave, spiralTime((c4(idx4)-1)*10+(1:10)')];
    forSave=[forSave,  fpcasTime((c5(idx5)-1)*10+(1:10)')];
    forSave=[forSave,  fistaTime((c6(idx6)-1)*10+(1:10)')];

    forSave=[forSave,    npgCost((c1(idx1)-1)*10+(1:10)')];
    forSave=[forSave,   npgcCost((c2(idx2)-1)*10+(1:10)')];
    forSave=[forSave,   npgsCost((c3(idx3)-1)*10+(1:10)')];
    forSave=[forSave, spiralCost((c4(idx4)-1)*10+(1:10)')];
    forSave=[forSave,  fpcasCost((c5(idx5)-1)*10+(1:10)')];
    forSave=[forSave,  fistaCost((c6(idx6)-1)*10+(1:10)')];

    forSave=[forSave,    npgRMSE((c1(idx1)-1)*10+(1:10)')];
    forSave=[forSave,   npgcRMSE((c2(idx2)-1)*10+(1:10)')];
    forSave=[forSave,   npgsRMSE((c3(idx3)-1)*10+(1:10)')];
    forSave=[forSave, spiralRMSE((c4(idx4)-1)*10+(1:10)')];
    forSave=[forSave,  fpcasRMSE((c5(idx5)-1)*10+(1:10)')];
    forSave=[forSave,  fistaRMSE((c6(idx6)-1)*10+(1:10)')];
    forSave=[forSave, snr(:)];
    forSave=[forSave, sparsaTime((c7(idx7)-1)*10+(1:10)')];
    forSave=[forSave, sparsaCost((c7(idx7)-1)*10+(1:10)')];
    forSave=[forSave, sparsaRMSE((c7(idx7)-1)*10+(1:10)')];
    forSave=[forSave, sparsnTime((c8(idx8)-1)*10+(1:10)')];
    forSave=[forSave, sparsnCost((c8(idx8)-1)*10+(1:10)')];
    forSave=[forSave, sparsnRMSE((c8(idx8)-1)*10+(1:10)')];
    save('varySNR.data','forSave','-ascii');

    as=1:5; idx=1:3:10;
    for mIdx=idx
        figure(900);
        semilogy(log10(u(mIdx))+as-3,    npgRMSE(mIdx,as),'r-*'); hold on;
        semilogy(log10(u(mIdx))+as-3,   npgcRMSE(mIdx,as),'r.-');
        semilogy(log10(u(mIdx))+as-3, sparsnRMSE(mIdx,as),'r-s');
        semilogy(log10(u(mIdx))+as-3, spiralRMSE(mIdx,as),'r-^');
        semilogy(log10(u(mIdx))+as-3,   npgsRMSE(mIdx,as),'k-s');
        semilogy(log10(u(mIdx))+as-3,  fpcasRMSE(mIdx,as),'g-o');
        semilogy(log10(u(mIdx))+as-3,  fistaRMSE(mIdx,as),'g-.');
        semilogy(log10(u(mIdx))+as-3, sparsaRMSE(mIdx,as),'g->');

        figure;
        semilogy(log10(u(mIdx))+as-3,    npgTime(mIdx,as),'r-*'); hold on;
        semilogy(log10(u(mIdx))+as-3,   npgcTime(mIdx,as),'r.-');
        semilogy(log10(u(mIdx))+as-3, sparsnTime(mIdx,as),'r-s');
        semilogy(log10(u(mIdx))+as-3, spiralTime(mIdx,as),'r-^');
        semilogy(log10(u(mIdx))+as-3,   npgsTime(mIdx,as),'k-s');
        semilogy(log10(u(mIdx))+as-3,  fpcasTime(mIdx,as),'g-o');
        semilogy(log10(u(mIdx))+as-3,  fistaTime(mIdx,as),'g-.');
        semilogy(log10(u(mIdx))+as-3, sparsaTime(mIdx,as),'g->');
        legend('npg','npgc','sparsn','spiral','npgs','fpcas','fista','sparas');
        title(sprintf('mIdx=%d',mIdx));
    end
    figure(900); legend('npg','npgc','sparsn','spiral','npgs','fpcas','fista','sparas');

    system(['mv varySNRTable.tex varySNR.data ' paperDir]);
end

if(any(runList==006))
    filename = [mfilename '_006.mat']; load(filename);

    K=5;
    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    fprintf('Poisson example\n');

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    spiralTime= mean(Cell.getField(spiral(:,:,1:K),'time'),3);
    fistaTime = mean(Cell.getField( fista(:,:,1:K),'time'),3);
    istTime   = mean(Cell.getField(   ist(:,:,1:K),'time'),3);

    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    spiralCost= mean(Cell.getField(spiral(:,:,1:K),'cost'),3);
    fistaCost = mean(Cell.getField( fista(:,:,1:K),'cost'),3);
    istCost   = mean(Cell.getField(   ist(:,:,1:K),'cost'),3);

    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);
    spiralRMSE= mean(Cell.getField(spiral(:,:,1:K),'reconerror'),3);
    fistaRMSE = mean(Cell.getField( fista(:,:,1:K),'RMSE'),3);
    istRMSE   = mean(Cell.getField(   ist(:,:,1:K),'RMSE'),3);

    [r,c1]=find(   npgRMSE==repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c3]=find( fistaRMSE==repmat(min( fistaRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c4]=find(spiralRMSE==repmat(min(spiralRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c5]=find(   istRMSE==repmat(min(   istRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find(  npgcRMSE==repmat(min(  npgcRMSE,[],2),1,5)); [r,idx6]=sort(r);
    disp([c1(idx1) ,c2(idx2) ,c3(idx3) ,c4(idx4) ,c5(idx5) ,c6(idx6)])
    idx1 = 3;    
    idx2 = 3;
    idx3 = 5;
    idx4 = 3;
    idx5 = 3;
    idx6 = 3;

    figure;
    semilogy(m,   npgRMSE(:,idx1),'r-*'); hold on;
    loglog(m,  npgsRMSE(:,idx2),'c-p');
    %loglog(m, fistaRMSE(:,idx3),'g-s');
    loglog(m,spiralRMSE(:,idx4),'k-^');
    %loglog(m,   istRMSE(:,idx5),'b-.');
    loglog(m,npgcRMSE(:,idx6),'k*-.');

    figure;
    plot(m,   npgTime(:,idx1),'r-*'); hold on;
    loglog(m,  npgsTime(:,idx2),'c-p');
    %loglog(m, fistaTime(:,idx3),'g-s');
    loglog(m,spiralTime(:,idx4),'k-^');
    %loglog(m,   istTime(:,idx5),'b-.');
    loglog(m,npgcTime(:,idx6),'k*-.');

    keyboard;

    forSave=[];
    forSave=[forSave,    npgTime(:,idx1)];
    forSave=[forSave,   npgsTime(:,idx2)];
    forSave=[forSave,  fistaTime(:,idx3)];
    forSave=[forSave, spiralTime(:,idx4)];
    forSave=[forSave,    istTime(:,idx5)];

    forSave=[forSave,    npgCost(:,idx1)];
    forSave=[forSave,   npgsCost(:,idx2)];
    forSave=[forSave,  fistaCost(:,idx3)];
    forSave=[forSave, spiralCost(:,idx4)];
    forSave=[forSave,    istCost(:,idx5)];

    forSave=[forSave,    npgRMSE(:,idx1)];
    forSave=[forSave,   npgsRMSE(:,idx2)];
    forSave=[forSave,  fistaRMSE(:,idx3)];
    forSave=[forSave, spiralRMSE(:,idx4)];
    forSave=[forSave,    istRMSE(:,idx5)];
    forSave=[forSave, m(:)];
    save('varyMeasurementPoisson.data','forSave','-ascii');

    mIdx=5; forSave=[]; t=0;
    t=t+1; temp=   npg{mIdx,idx1,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=  npgs{mIdx,idx2,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp= fista{mIdx,idx3,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    t=t+1; temp=   ist{mIdx,idx5,1}.stepSize(:); forSave(1:length(temp),t)=temp;
    save('stepSize.data','forSave','-ascii');

    forSave=[]; t=0;
    out=   ist{mIdx,idx5,1};
    t=t+1; forSave(1:length(out.cost      ),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE      ),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time      ),t)=out.time;
    out=   npg{mIdx,idx1,1};
    t=t+1; forSave(1:length(out.cost      ),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE      ),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time      ),t)=out.time;
    out=spiral{mIdx,idx4,1};
    t=t+1; forSave(1:length(out.cost      ),t)=out.cost;
    t=t+1; forSave(1:length(out.reconerror),t)=out.reconerror;
    t=t+1; forSave(1:length(out.time      ),t)=out.time;

    keyboard
    mincost=reshape(forSave(:,[1,4,7]),[],1); 
    mincost=min(mincost(mincost~=0));
    idx=(forSave(:,1)~=0); forSave(idx,1)=(forSave(idx,1)-mincost);
    idx=(forSave(:,4)~=0); forSave(idx,4)=(forSave(idx,4)-mincost);
    idx=(forSave(:,7)~=0); forSave(idx,7)=(forSave(idx,7)-mincost);


    save('cost_itr.data','forSave','-ascii');
    system(['mv varyMeasurementPoisson.data stepSize.data cost_itr.data ' paperDir]);
end

if(any(runList==007))
    filename = [mfilename '_007.mat']; load(filename);

    m=[ 200, 300, 400, 500, 600, 700, 800]; % should go from 200
    u=[1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2];
    idx=2:2:7;
    K = 5;

    npgTime   = mean(Cell.getField(   npg(:,:,1:K),'time'),3);
    npgCost   = mean(Cell.getField(   npg(:,:,1:K),'cost'),3);
    npgRMSE   = mean(Cell.getField(   npg(:,:,1:K),'RMSE'),3);

    npgcTime  = mean(Cell.getField(  npgc(:,:,1:K),'time'),3);
    npgcCost  = mean(Cell.getField(  npgc(:,:,1:K),'cost'),3);
    npgcRMSE  = mean(Cell.getField(  npgc(:,:,1:K),'RMSE'),3);

    npgsTime  = mean(Cell.getField(  npgs(:,:,1:K),'time'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,:,1:K),'cost'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,:,1:K),'RMSE'),3);

    npgscTime = mean(Cell.getField( npgsc(:,:,1:K),'time'),3);
    npgscCost = mean(Cell.getField( npgsc(:,:,1:K),'cost'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,:,1:K),'RMSE'),3);

    npgsclTime= mean(Cell.getField(npgscl(:,:,1:K),'time'),3);
    npgsclCost= mean(Cell.getField(npgscl(:,:,1:K),'cost'),3);
    npgsclRMSE= mean(Cell.getField(npgscl(:,:,1:K),'RMSE'),3);

    npgclTime = mean(Cell.getField( npgcl(:,:,1:K),'time'),3);
    npgclCost = mean(Cell.getField( npgcl(:,:,1:K),'cost'),3);
    npgclRMSE = mean(Cell.getField( npgcl(:,:,1:K),'RMSE'),3);

    for i=1:length(m)
        for j=1:size(npgTime,2)
            for k=1:K
                a = log10(u(i))+(j-3);
                temp = abs(gnet{i,k}.a-a);
                idx = find(temp==min(temp));
                gnetTime(i,j,k)=gnet{i,k}.time;
                gnetCost(i,j,k)=gnet{i,k}.cost(idx);
                gnetRMSE(i,j,k)=gnet{i,k}.RMSE(idx);
            end
            gnetAlpha{i,j}=gnet{i,1}.alpha(:,idx)
        end
    end
    gnetTime = mean(gnetTime,3);
    gnetCost = mean(gnetCost,3);
    gnetRMSE = mean(gnetRMSE,3);

    [r,c1]=find(    npgRMSE == repmat(min(   npgRMSE,[],2),1,5)); [r,idx1]=sort(r);
    [r,c2]=find(   npgcRMSE == repmat(min(  npgcRMSE,[],2),1,5)); [r,idx2]=sort(r);
    [r,c3]=find(   npgsRMSE == repmat(min(  npgsRMSE,[],2),1,5)); [r,idx3]=sort(r);
    [r,c4]=find(  npgscRMSE == repmat(min( npgscRMSE,[],2),1,5)); [r,idx4]=sort(r);
    [r,c5]=find( npgsclRMSE == repmat(min(npgsclRMSE,[],2),1,5)); [r,idx5]=sort(r);
    [r,c6]=find(  npgclRMSE == repmat(min( npgclRMSE,[],2),1,5)); [r,idx6]=sort(r);

    disp([c1(idx1), c2(idx2), c3(idx3), c4(idx4), c5(idx5), c6(idx6)]);

    figure;
    semilogy(m,   npgRMSE((c1(idx1)-1)*7+(1:7)'),'r-*'); hold on;
    semilogy(m,  npgcRMSE((c2(idx2)-1)*7+(1:7)'),'r-p');
    semilogy(m,  npgsRMSE((c3(idx3)-1)*7+(1:7)'),'k-s');
    semilogy(m, npgscRMSE((c4(idx4)-1)*7+(1:7)'),'k-<');
    semilogy(m, npgclRMSE((c5(idx5)-1)*7+(1:7)'),'r:p');
    semilogy(m,npgsclRMSE((c6(idx6)-1)*7+(1:7)'),'k:<');
    semilogy(m,  gnetRMSE((c4(idx4)-1)*7+(1:7)'),'c-o');
    legend('npg','npgc','npgs','npgsc','npgcl','npgscl','gnet');
    figure;
    semilogy(m,   npgTime((c1(idx1)-1)*7+(1:7)'),'r-*'); hold on;
    semilogy(m,  npgcTime((c2(idx2)-1)*7+(1:7)'),'r-p');
    semilogy(m,  npgsTime((c3(idx3)-1)*7+(1:7)'),'k-s');
    semilogy(m, npgscTime((c4(idx4)-1)*7+(1:7)'),'k-<');
    semilogy(m, npgclTime((c5(idx5)-1)*7+(1:7)'),'r:p');
    semilogy(m,npgsclTime((c6(idx6)-1)*7+(1:7)'),'k:<');
    semilogy(m,  gnetTime((c4(idx4)-1)*7+(1:7)'),'c-o');
    legend('npg','npgc','npgs','npgsc','npgcl','npgscl','gnet');

    keyboard

    forSave=[];
    forSave=[forSave,    npgTime((c1(idx1)-1)*7+(1:7)')];
    forSave=[forSave,    npgCost((c1(idx1)-1)*7+(1:7)')];
    forSave=[forSave,    npgRMSE((c1(idx1)-1)*7+(1:7)')];
    
    forSave=[forSave,   npgcTime((c2(idx2)-1)*7+(1:7)')];
    forSave=[forSave,   npgcCost((c2(idx2)-1)*7+(1:7)')];
    forSave=[forSave,   npgcRMSE((c2(idx2)-1)*7+(1:7)')];

    forSave=[forSave,   npgsTime((c3(idx3)-1)*7+(1:7)')];
    forSave=[forSave,   npgsCost((c3(idx3)-1)*7+(1:7)')];
    forSave=[forSave,   npgsRMSE((c3(idx3)-1)*7+(1:7)')];

    forSave=[forSave,  npgscTime((c4(idx4)-1)*7+(1:7)')];
    forSave=[forSave,  npgscCost((c4(idx4)-1)*7+(1:7)')];
    forSave=[forSave,  npgscRMSE((c4(idx4)-1)*7+(1:7)')];

    forSave=[forSave,  npgclTime((c5(idx5)-1)*7+(1:7)')];
    forSave=[forSave,  npgclCost((c5(idx5)-1)*7+(1:7)')];
    forSave=[forSave,  npgclRMSE((c5(idx5)-1)*7+(1:7)')];

    forSave=[forSave, npgsclTime((c6(idx6)-1)*7+(1:7)')];
    forSave=[forSave, npgsclCost((c6(idx6)-1)*7+(1:7)')];
    forSave=[forSave, npgsclRMSE((c6(idx6)-1)*7+(1:7)')];

    forSave=[forSave,   gnetTime((c4(idx4)-1)*7+(1:7)')];
    forSave=[forSave,   gnetCost((c4(idx4)-1)*7+(1:7)')];
    forSave=[forSave,   gnetRMSE((c4(idx4)-1)*7+(1:7)')];

    forSave=[forSave, m(:)];
    save('varyMeasurementSkylineLogPoisson.data','forSave','-ascii');

    forSave=m(:);
    forSave=[forSave,    npgTime((c1(idx1)-1)*7+(1:7)')];
    forSave=[forSave,   npgcTime((c2(idx2)-1)*7+(1:7)')];
    forSave=[forSave,   npgsTime((c3(idx3)-1)*7+(1:7)')];
    forSave=[forSave,  npgscTime((c4(idx4)-1)*7+(1:7)')];
    forSave=[forSave,  npgclTime((c5(idx5)-1)*7+(1:7)')];
    forSave=[forSave, npgsclTime((c6(idx6)-1)*7+(1:7)')];
    forSave=[forSave,   gnetTime((c4(idx4)-1)*7+(1:7)')];
    save('varyMeasurementTimeSkylineLogPoisson.data','forSave','-ascii');

    return;

    f=fopen('selectedTime.data','w');
    for mIdx=1:length(m)
        fprintf(f,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t%s\t%s\n',...
                npgTime(mIdx,uNonneg(mIdx)), ...
               npgcTime(mIdx,uNonneg(mIdx)), ...
             spiralTime(mIdx,uNonneg(mIdx)), ...
             sparsnTime(mIdx,uNonneg(mIdx)), ...
               npgsTime(mIdx,uNeg   (mIdx)), ...
              fpcasTime(mIdx,uNeg   (mIdx)), ...
              fistaTime(mIdx,uNeg   (mIdx)), ...
             sparsaTime(mIdx,uNeg   (mIdx)), ...
            m(mIdx),num2str(log10(u((mIdx)))+uNonneg(mIdx)-3), num2str(log10(u(mIdx))+uNeg(mIdx)-3));
    end
    fclose(f);

    temp = 4;
    signal=npg{1}.opt.trueAlpha;
    signal=[signal,    npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,   npgc{gEle((c2(idx2)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,   npgs{gEle((c3(idx3)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, spiral{gEle((c4(idx4)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  fpcas{gEle((c5(idx5)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal,  fista{gEle((c6(idx6)-1)*9+(1:9)',temp)}.alpha];
    signal=[signal, sparsa{gEle((c7(idx7)-1)*9+(1:9)',temp)}.alpha];
    save('skyline.data','signal','-ascii');
    figure; plot(signal(:,2)); hold on; plot(signal(:,1),'r'); title('NPG');
    figure; plot(signal(:,4)); hold on; plot(signal(:,1),'r'); title('NPGs');
    figure; plot(signal(:,6)); hold on; plot(signal(:,1),'r'); title('FPCas');
    fprintf('npgRec RMSE: %g%%\n',npg{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100);
    fprintf('npgsRec RMSE: %g%%\n',npgs{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100);
    fprintf('fistaRec RMSE: %g%%\n',fista{gEle((c1(idx1)-1)*9+(1:9)',temp)}.RMSE(end)*100);

end

if(any(runList==008))
    filename = [mfilename '_008.mat']; load(filename);
    i=3; t=0;
    prjFull = [60, 80, 100, 120, 180, 360];
    forSave=[];

    out=   npg{i};
    fprintf('NPGwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'NPGwLin.png','png');
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    out=spiral{i};
    fprintf('SPIRALwLin: %g\n',out.reconerror(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'SPIRALwLin.png','png');
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.reconerror),t)=out.reconerror;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    out=fbp{i};
    fprintf('FBPwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,spiral{i}.opt.mask);
    imwrite(img/max(img(:)),'FBPwLin.png','png');

    out=fpcas{i};
    fprintf('FPCASwLin: %g\n',out.RMSE(end));
    img=showImgMask(out.alpha,out.opt.mask);
    imwrite(img/max(img(:)),'FPCASwLin.png','png');

    perform=[]; a=1;
    perform=[perform, reshape(Cell.getField(fbp      (:,1),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(fbp      (:,1),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(npg      (:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(npg      (:,a),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(spiral   (:,1),'reconerror'),[],1)];
    perform=[perform, reshape(Cell.getField(spiral   (:,1),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(fpcas    (:,1),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(fpcas    (:,1),'cpu' ),[],1)];
    perform=[perform, reshape(Cell.getField(npgc     (:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(npgc     (:,a),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(sparsa   (:,1),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(sparsa   (:,1),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(npgs     (:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(npgs     (:,a),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(npgsc    (:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(npgsc    (:,a),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(sparsn   (:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(sparsn   (:,a),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(npgc_nads(:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(npgc_nads(:,a),'time'),[],1)];
    perform=[perform, reshape(Cell.getField(npg_nads (:,a),'RMSE'),[],1)];
    perform=[perform, reshape(Cell.getField(npg_nads (:,a),'time'),[],1)];

    figure;
    style1 = {'b-','r-','b:','r:','b--','r--','-','*','.','+','d'};
    style2 = {'b-*','r-*','b:o','r:o','b--s','r--s','c^-','k-.p','c:>','g-.p','k.-'};
    for i=1:9
        semilogy(prjFull/2,perform(:,i*2-1),style2{i}); hold on;
    end
    legend('fbp','npg','spiral','fpcas','npgc','sparsa','npgs','npgsc','sparsn','npgc_nads','npg_nads');
    figure; 
    for i=1:9
        semilogy(prjFull/2,perform(:,i*2),style2{i}); hold on;
    end
    legend('fbp','npg','spiral','fpcas','npgc','sparsa','npgs','npgsc','sparsn','npgc_nads','npg_nads');

    keyboard

    perform=[perform,prjFull(:)];
    save('varyPrj.data','perform','-ascii');
    save('xray_itr.data','forSave','-ascii');
    !cp varyPrj.data xray_itr.data *wLin.png ~/research/myPaper/asilomar2014/
end

if(any(runList==009))
    filename = [mfilename '_009.mat']; load(filename);

    prjFull = [60, 80, 100, 120, 180, 360]; j=1;
    fprintf('Poisson Log link example with glass beads\n');

    k=1;
    npgTime   = showResult(   npg(:,:,k),2,'time');
    npgsTime  = showResult(  npgs(:,:,k),2,'time');
    % fbpTime   = showResult(   fbp(:,:,k),2,'time');

    npgCost   = showResult(   npg(:,:,k),2,'cost');
    npgsCost  = showResult(  npgs(:,:,k),2,'cost');
    % fbpCost   = showResult(   fbp(:,:,k),2,'cost');

    npgRMSE   = showResult(   npg(:,:,k),2,'RMSE');
    npgsRMSE  = showResult(  npgs(:,:,k),2,'RMSE');
    fbpRMSE   = showResult(   fbp(:,:,k),2,'RMSE');

    % [r,c1]=find(   npgRMSE==repmat(min(   npgRMSE,[],2),1,3)); [r,idx1]=sort(r);
    % [r,c2]=find(  npgsRMSE==repmat(min(  npgsRMSE,[],2),1,3)); [r,idx2]=sort(r);
    % [r,c3]=find(   fbpRMSE==repmat(min(   fbpRMSE,[],2),1,3)); [r,idx3]=sort(r);
    % [c1(idx1) ,c2(idx2) ,c3(idx3)]
    idx1 = 2;    
    idx2 = 2;
    idx3 = 2;

    figure;
    semilogy(prjFull,   npgRMSE(:,idx1),'r-*'); hold on;
    loglog(prjFull,  npgsRMSE(:,idx2),'c-p');
    loglog(prjFull, fbpRMSE(:,idx3),'g-s');

    figure;
    plot(prjFull,   npgTime(:,idx1),'r-*'); hold on;
    loglog(prjFull,  npgsTime(:,idx2),'c-p');
    %loglog(prjFull,   fbpTime(:,idx3),'g-s');

    forSave=[];
    forSave=[forSave,    npgRMSE(:,idx1)];
    forSave=[forSave,   npgsRMSE(:,idx2)];
    forSave=[forSave,    fbpRMSE(:,idx3)];

    forSave=[forSave,    npgTime(:,idx1)];
    forSave=[forSave,   npgsTime(:,idx2)];
    % forSave=[forSave,    fbpTime(:,idx3)];

    % forSave=[forSave,    npgCost(:,idx1)];
    % forSave=[forSave,   npgsCost(:,idx2)];
    % forSave=[forSave,    fbpCost(:,idx3)];

    forSave=[forSave, prjFull(:)];
    save('varyPrjGlassBead.data','forSave','-ascii');

    idx=4;
    img=showImgMask(npg {idx,idx1}.alpha,npg{idx,idx1}.opt.mask);
    imwrite(img/max(img(:)),'NPGgb.png','png');
    img=showImgMask(npgs{idx,idx2}.alpha,npg{idx,idx1}.opt.mask);
    imwrite(img/max(img(:)),'NPGSgb.png','png');
    img=showImgMask(fbp {idx,idx3}.alpha,npg{idx,idx1}.opt.mask);
    imwrite(img/max(img(:)),'FBPgb.png','png');

    disp([npg{idx,idx1}.RMSE(end), npgs{idx,idx2}.RMSE(end), fbp{idx,idx3}.RMSE]);

    system(['mv varyPrjGlassBead.data NPGgb.png NPGSgb.png FBPgb.png ' paperDir]);
    keyboard
end
end

function e=gEle(x,i); e=x(i); end

