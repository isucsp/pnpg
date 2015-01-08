function runAsilomar2014(runList)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Beam Hardening correction of CT Imaging via Mass attenuation 
%                        coefficient discretizati
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration

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

end


