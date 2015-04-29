function yangEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polychromatic Sparse Image Reconstruction and Mass Attenuation Spectrum 
%            Estimation via B-Spline Basis Function Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:  Changed to class oriented for easy configuration

% dis discretized, polychromatic model
% dis, compare the NPG and NCG_PR for both continuation and 
% non-continuation
%
% compare the effect forcing center model of spectrum 
%
% A few thing to tune: continuation, centerb, ActiveSet VS FISTA_Simplex
% maxIeSteps

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt'); filename = [mfilename '.mat'];

        opt.beamharden=true; opt.spectBasis='b1'; opt.E=30;
        opt.estIe=true;

        prjFull = [60, 80, 100, 120, 180, 360];
        for i=length(prjFull)
            opt.prjFull = prjFull(i); opt.prjNum = opt.prjFull;

            [y,Phi,Phit,Psi,Psit,opt,FBP]=loadYang(opt);
            opt.maxItr=4e3; opt.thresh=1e-6;

            initSig = maskFunc(FBP(y),opt.mask~=0);

            j=1;
            fprintf('%s, i=%d, j=%d\n','Filtered Backprojection',i,j);
            fbp{i}.img=FBP(y);
            fbp{i}.alpha=fbp{i}.img(opt.mask~=0);
            fbp{i}.RMSE=1-(innerProd(fbp{i}.alpha,opt.trueAlpha)^2)/sqrNorm(opt.trueAlpha)/sqrNorm(fbp{i}.alpha);
            fprintf('fbp RMSE=%g\n',fbp{i}.RMSE);

            % unknown ι(κ), NPG-AS
            for j=[2:4]
                fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
                u  =  10.^[-5  -5   -5   -5   -5   -5];
                opt.u=u(i)*10^(j-3);
                maxItr=2000;
                npg2_b1{i,j}=BHC.NPG2(Phi,Phit,Psi,Psit,y,initSig,opt);

%               fpcas {i,j}=Wrapper.FPCas(Phi,Phit,Psi,Psit,y,initSig,opt);
                save(filename);
            end

            continue;

            % known ι(κ), 

            aArray=[-10:-2];
            for j=5:7
                fprintf('%s, i=%d, j=%d\n','FPCAS',i,j);

                opt.u=10^aArray(j)*max(abs(At(y)));
                FPC_AS
            end
        end

    case 'plot'
        conf=ConfigCT('castSim','CircleMask','gpuPrj');
        opt.prjFull = 360; opt.prjNum = opt.prjFull;
        y = Phi(opt.trueAlpha); % equivalent to linear projection
        forSave=[conf.y y];
        polyy=linspace(min(conf.y),max(conf.y),1000);
        [y1,idx]=sort(conf.y); y2=y(idx);
        [y1,ia,ic]=unique(y1); y2=y2(ia);
        lineary=interp1(y1,y2,polyy,'spline');
        forsave=[polyy(:); lineary(:)];
        save('linearization.data','forSave','-ascii');

        return;

        load([mfilename '_021.mat']);
        fprintf('for non skiped Ie\n');
        t=1; noSkiped(:,t)=1:2000;
        for i=[1,3]
            t=t+1; noSkiped(1:length(out{1,i}.cost),t)=out{1,i}.cost;
            t=t+1; noSkiped(1:length(out{1,i}.RMSE),t)=out{1,i}.RMSE;
            t=t+1; noSkiped(1:length(out{1,i}.time),t)=out{1,i}.time;
        end
        mincost=reshape(noSkiped(:,[2,5]),[],1); mincost=min(mincost(mincost>0));
        noSkiped(:,[2,5])=noSkiped(:,[2,5])-mincost;
        save('costRmseTime.data','noSkiped','-ascii');
        return;
        figure; subplot(2,1,1);
        semilogy(out{1,1}.cost-mincost,'r'); hold on;
        semilogy(out{1,2}.cost-mincost,'g');
        semilogy(out{1,3}.cost-mincost,'b'); subplot(2,1,2);
        semilogy(out{1,1}.RMSE,'r'); hold on;
        semilogy(out{1,2}.RMSE,'g');
        semilogy(out{1,3}.RMSE,'b');
        figure; subplot(2,1,1);
        semilogy(out{1,1}.time, out{1,1}.cost-mincost,'r'); hold on;
        semilogy(out{1,2}.time, out{1,2}.cost-mincost,'g');
        semilogy(out{1,3}.time, out{1,3}.cost-mincost,'b'); subplot(2,1,2);
        semilogy(out{1,1}.time, out{1,1}.RMSE,'r'); hold on;
        semilogy(out{1,2}.time, out{1,2}.RMSE,'g');
        semilogy(out{1,3}.time, out{1,3}.RMSE,'b');

        return;
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

        polymodel = Spline(out1.opt.spectBasis,out1.kappa);
        polymodel.setPlot(out1.opt.trueKappa,out1.opt.trueIota,out1.opt.epsilon);
        [tK,tU,kappa,Ie]=polymodel.plotSpectrum(out1.Ie);
        idx=find(tU==max(tU));
        q=tK(idx)/1; q=1/1.15
        tK=tK/q; tU=tU*q;
        i=1; forSave=[tK, tU];

        i=i+2; forSave(1:length(out1.kappa),i:i+1)=[out1.kappa,out1.Ie];
        figure; semilogx(tK,tU,'r'); hold on; loglog(out1.kappa,out1.Ie,'g-*');

        i=i+2; forSave(1:length(out2.kappa),i:i+1)=[out2.kappa,out2.Ie];
        loglog(out2.kappa,out2.Ie,'b-<');

        i=i+2; forSave(1:length(out3.kappa),i:i+1)=[out3.kappa,out3.Ie];
        loglog(out3.kappa,out3.Ie,'c-s');
        save('effectiveCenterB.data','forSave','-ascii');
        
        paperDir = './';
        %system(['mv effectiveCenterB.data ' paperDir]);
end
end


