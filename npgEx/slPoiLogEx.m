function slPoiLogEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration
%
%                    Skyline log link Poisson Example
%    vary the number of measurements and inital intensity constatn I_0

if(~exist('op','var')) op='plot'; end

switch lower(op)
    case {lower('I016'), lower('I012'), lower('I009')}
        if(strcmpi(op,'I016'))
            filename = [mfilename '_007.mat'];
            if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end;
            clear('opt'); opt.B=2^16; aa=-3:-0.1:-6;
        elseif(strcmpi(op,'I012'))
            filename = [mfilename '_027.mat'];
            if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end;
            clear('opt'); opt.B=2^12; aa=-1:-0.2:-5;
        elseif(strcmpi(op,'I009'))
            filename = [mfilename '_037.mat'];
            if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end;
            clear('opt'); opt.B=2^9; aa=-0.2:-0.2:-4;
        end

        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e5; opt.thresh=1e-6;
        m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200
        a=-4.4:-0.2:-4.8;
        for k=1:10
            for i=1:length(m)
                opt.m=m(i); opt.noiseType='poissonLogLink'; opt.matrixType='conv';
                [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
                initSig=-Phit(invEAAt*log(max(y,1)/max(y)))*0;
                fprintf('i=%d, k=%d, min=%d, max=%d\n',i,k,min(y), max(y));

                %% fit by approximated Gaussian model for known I0
                % u_max=pNorm(Psit(Phit(y.*log(opt.I0./max(y,1)))),inf);
                %% fit by Poisson model with known I0
                % u_max=pNorm(Psit(Phit(y-opt.I0)),inf);
                %% fit by the unkown I0 log link Poisson model
                u_max=pNorm(Psit(Phit(y-mean(y))),inf);

                opt.noiseType='poissonLogLink';
                temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
                xx=opt.thresh; opt.thresh=1e-12;
                gnet    {i,k}=Wrapper.glmnet(Utils.getMat(Phi,length(initSig)),Utils.getMat(Psi,length(Psit(initSig))),y,initSig,opt);
                opt.thresh=xx; clear('xx');
                npgFull {i,k}=Wrapper.NPG  (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgsFull{i,k}=Wrapper.NPGs (Phi,Phit,Psi,Psit,y,initSig,opt);
                opt=temp;

                opt.noiseType='poissonLogLink0';
                temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
                npgFull_knownI0 {i,k}=Wrapper.NPG (Phi,Phit,Psi,Psit,y,initSig,opt);
                npgsFull_knownI0{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt);
                xx=opt.thresh; opt.thresh=1e-12;
                gnet0   {i,k}=Wrapper.glmnet(Utils.getMat(Phi,length(initSig)),Utils.getMat(Psi,length(Psit(initSig))),y,initSig,opt);
                opt.thresh=xx; clear('xx');
                opt=temp;

%               opt.noiseType='gaussian';
%               %% fit by approximated Gaussian model for known I0
%               yy=log(opt.I0./max(y,1)); yy=yy.*sqrt(y);
%               wPhi=@(xxx) sqrt(y).*Phi(xxx); wPhit=@(xxx) Phit(sqrt(y).*xxx);

%               temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
%               npglwFull {i,k}=Wrapper.NPG(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
%               npgslwFull{i,k}=Wrapper.NPGs(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
%               opt=temp;

                save(filename);
            end
        end

    case lower('ind')

        filename = [mfilename '_017.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt');
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e5; opt.thresh=1e-6;
        m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200
        a=-4.4:-0.2:-4.8;
        for k=1:5
            for i=1:length(m)
                opt.m=m(i); opt.noiseType='poissonLogLink'; opt.matrixType='conv';
                [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
                initSig=-Phit(invEAAt*log(max(y,1)/max(y)))*0;
                fprintf('i=%d, k=%d, min=%d, max=%d\n',i,k,min(y), max(y));

                %% fit by approximated Gaussian model for known I0
                % u_max=pNorm(Psit(Phit(y.*log(opt.I0./max(y,1)))),inf);
                yy=log(opt.I0./max(y,1)); yy=yy.*sqrt(y);
                wPhi=@(xxx) sqrt(y).*Phi(xxx); wPhit=@(xxx) Phit(sqrt(y).*xxx);
                %% fit by Poisson model with known I0
                % u_max=pNorm(Psit(Phit(y-opt.I0)),inf);
                %% fit by the unkown I0 log link Poisson model
                u_max=pNorm(Psit(Phit(y-mean(y))),inf);

                for j=1:length(a)
                    opt.u=10^a(j)*u_max;

                    opt.noiseType='poissonLogLink';
                    npg {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npgs{i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);

                    opt.noiseType='poissonLogLink0';
                    npg0 {i,j,k}=Wrapper.NPGc (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npgs0{i,j,k}=Wrapper.NPGsc(Phi,Phit,Psi,Psit,y,initSig,opt);

                    opt.noiseType='gaussian';
                    npglw {i,j,k}=Wrapper.NPGc(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
                    npgslw{i,j,k}=Wrapper.NPGsc(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
                end
                save(filename);
            end
        end

    case 'plot1'
        filename = [mfilename '_007.mat']; load(filename);

        k=1;
        m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200

        for i=1:length(m)
            npgRMSE(i,k) = min(npgFull{i,k}.contRMSE);
            npgsRMSE(i,k) = min(npgsFull{i,k}.contRMSE);
            gnetRMSE(i,k) = min(gnet{i,k}.RMSE(:));
            gnet0RMSE(i,k) = min(gnet0{i,k}.RMSE(:));
            npg0RMSE(i,k) = min(npgFull_knownI0{i,k}.contRMSE);
            npgs0RMSE(i,k) = min(npgsFull_knownI0{i,k}.contRMSE);
            npglwRMSE(i,k) = min(npglwFull{i,k}.contRMSE);
            npgslwRMSE(i,k) = min(npgslwFull{i,k}.contRMSE);

            npgRec(i,k)=         aa(min(find(          npgFull{i,k}.contRMSE==   npgRMSE(i,k))));
            npgsRec(i,k)=         aa(min(find(         npgsFull{i,k}.contRMSE==  npgsRMSE(i,k))));
            npg0Rec(i,k)=         aa(min(find(  npgFull_knownI0{i,k}.contRMSE==  npg0RMSE(i,k))));
            npgs0Rec(i,k)=         aa(min(find( npgsFull_knownI0{i,k}.contRMSE== npgs0RMSE(i,k))));
            npglwRec(i,k)=         aa(min(find(        npglwFull{i,k}.contRMSE== npglwRMSE(i,k))));
            npgslwRec(i,k)=         aa(min(find(       npgslwFull{i,k}.contRMSE==npgslwRMSE(i,k))));
            glmnetRec(i,k)=gnet{i,k}.a(min(find(             gnet{i,k}.RMSE    ==  gnetRMSE(i,k))));
            glmnet0Rec(i,k)=gnet0{i,k}.a(min(find(             gnet0{i,k}.RMSE    ==  gnet0RMSE(i,k))));
        end

        figure;
        semilogy(m,   npgRMSE(:,k),'r-*'); hold on;
        plot(    m,  npgsRMSE(:,k),'b-*');
        plot(    m,  npg0RMSE(:,k),'r.-');
        plot(    m, npgs0RMSE(:,k),'b.-');
        plot(    m,  gnetRMSE(:,k),'ch-');
        plot(    m, gnet0RMSE(:,k),'kh-');
        % plot(    m,  npglRMSE,'rp-');
        % plot(    m, npgslRMSE,'bp-');
        legend('npg','npgs','npg0','npgs0','glmnet','glmnet0');

        forSave=[ m(:), npgRMSE(:,k), ...
            npgsRMSE(:,k), ...
            npg0RMSE(:,k), ...
            npgs0RMSE(:,k), ...
            npglwRMSE(:,k), ...
            npgslwRMSE(:,k), ...
            gnetRMSE(:,k), ...
            gnet0RMSE(:,k)];

        disp([mean(    npgRec,2), ...
            mean(   npg0Rec,2), ...
            mean(  npglwRec,2), ...
            mean(  npgs0Rec,2), ...
            mean( npgslwRec,2), ...
            mean(   npgsRec,2), ...
            mean( glmnetRec,2), ...
            mean(glmnet0Rec,2)]);

        filename = [mfilename '_027.mat']; load(filename);
        for i=1:length(m)
            npgRMSE(i,k) = min(npgFull{i,k}.contRMSE);
            npgsRMSE(i,k) = min(npgsFull{i,k}.contRMSE);
            gnetRMSE(i,k) = min(gnet{i,k}.RMSE(:));
            npg0RMSE(i,k) = min(npgFull_knownI0{i,k}.contRMSE);
            npgs0RMSE(i,k) = min(npgsFull_knownI0{i,k}.contRMSE);
            npglwRMSE(i,k) = min(npglwFull{i,k}.contRMSE);
            npgslwRMSE(i,k) = min(npgslwFull{i,k}.contRMSE);

            npgRec(i,k)=         aa(min(find(          npgFull{i,k}.contRMSE==   npgRMSE(i,k))));
            npgsRec(i,k)=         aa(min(find(         npgsFull{i,k}.contRMSE==  npgsRMSE(i,k))));
            npg0Rec(i,k)=         aa(min(find(  npgFull_knownI0{i,k}.contRMSE==  npg0RMSE(i,k))));
            npgs0Rec(i,k)=         aa(min(find( npgsFull_knownI0{i,k}.contRMSE== npgs0RMSE(i,k))));
            npglwRec(i,k)=         aa(min(find(        npglwFull{i,k}.contRMSE== npglwRMSE(i,k))));
            npgslwRec(i,k)=         aa(min(find(       npgslwFull{i,k}.contRMSE==npgslwRMSE(i,k))));
            glmnetRec(i,k)=gnet{i,k}.a(min(find(             gnet{i,k}.RMSE    ==  gnetRMSE(i,k))));
        end
        figure;
        semilogy(m,   npgRMSE(:,k),'r-*'); hold on;
        plot(    m,  npgsRMSE(:,k),'b-*');
        plot(    m,  npg0RMSE(:,k),'r.-');
        plot(    m, npgs0RMSE(:,k),'b.-');
        plot(    m,  gnetRMSE(:,k),'ch-');
        % plot(    m,  npglRMSE,'rp-');
        % plot(    m, npgslRMSE,'bp-');
        legend('npg','npgs','npg0','npgs0','glmnet');
        disp([mean(   npgRec,2), ...
            mean(  npg0Rec,2), ...
            mean( npglwRec,2), ...
            mean( npgs0Rec,2), ...
            mean(npgslwRec,2), ...
            mean(  npgsRec,2), ...
            mean(glmnetRec,2)]);

        filename = [mfilename '_037.mat']; load(filename);
        for i=1:length(m)
            npgRMSE(i,k) = min(npgFull{i,k}.contRMSE);
            npgsRMSE(i,k) = min(npgsFull{i,k}.contRMSE);
            gnetRMSE(i,k) = min(gnet{i,k}.RMSE(:));
            gnet0RMSE(i,k) = min(gnet0{i,k}.RMSE(:));
            npg0RMSE(i,k) = min(npgFull_knownI0{i,k}.contRMSE);
            npgs0RMSE(i,k) = min(npgsFull_knownI0{i,k}.contRMSE);
            npglwRMSE(i,k) = min(npglwFull{i,k}.contRMSE);
            npgslwRMSE(i,k) = min(npgslwFull{i,k}.contRMSE);

            npgRec(i,k)=         aa(min(find(          npgFull{i,k}.contRMSE==   npgRMSE(i,k))));
            npgsRec(i,k)=         aa(min(find(         npgsFull{i,k}.contRMSE==  npgsRMSE(i,k))));
            npg0Rec(i,k)=         aa(min(find(  npgFull_knownI0{i,k}.contRMSE==  npg0RMSE(i,k))));
            npgs0Rec(i,k)=         aa(min(find( npgsFull_knownI0{i,k}.contRMSE== npgs0RMSE(i,k))));
            npglwRec(i,k)=         aa(min(find(        npglwFull{i,k}.contRMSE== npglwRMSE(i,k))));
            npgslwRec(i,k)=         aa(min(find(       npgslwFull{i,k}.contRMSE==npgslwRMSE(i,k))));
            glmnetRec(i,k)=gnet{i,k}.a(min(find(             gnet{i,k}.RMSE    ==  gnetRMSE(i,k))));
            glmnet0Rec(i,k)=gnet0{i,k}.a(min(find(             gnet0{i,k}.RMSE    ==  gnet0RMSE(i,k))));
        end
        figure;
        semilogy(m,   npgRMSE(:,k),'r-*'); hold on;
        plot(    m,  npgsRMSE(:,k),'b-*');
        plot(    m,  npg0RMSE(:,k),'r.-');
        plot(    m, npgs0RMSE(:,k),'b.-');
        plot(    m,  gnetRMSE(:,k),'ch-');
        plot(    m, gnet0RMSE(:,k),'kh-');
        % plot(    m,  npglRMSE,'rp-');
        % plot(    m, npgslRMSE,'bp-');
        legend('npg','npgs','npg0','npgs0','glmnet','glmnet0');
        disp([mean(    npgRec,2), ...
            mean(   npg0Rec,2), ...
            mean(  npglwRec,2), ...
            mean(  npgs0Rec,2), ...
            mean( npgslwRec,2), ...
            mean(   npgsRec,2), ...
            mean( glmnetRec,2), ...
            mean(glmnet0Rec,2)]);

        forSave=[ forSave, npgRMSE(:,k), ...
            npgsRMSE(:,k), ...
            npg0RMSE(:,k), ...
            npgs0RMSE(:,k), ...
            npglwRMSE(:,k), ...
            npgslwRMSE(:,k), ...
            gnetRMSE(:,k),...
            gnet0RMSE(:,k)];
        save('rmse_I0N.data','forSave','-ascii');
        (1-sort(forSave(:,5)./forSave(:,3) )')*100
        (1-sort(forSave(:,13)./forSave(:,11))')*100

    case 'plot2'

        filename = [mfilename '_017.mat']; load(filename);
        ttt=load([mfilename '_007.mat']); gnet=ttt.gnet; clear('ttt');

        m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200

        npgRMSE = mean(Cell.getField(    npg,'RMSE'),3);
        npgsRMSE = mean(Cell.getField(   npgs,'RMSE'),3);
        npg0RMSE = mean(Cell.getField(   npg0,'RMSE'),3);
        npgs0RMSE = mean(Cell.getField(  npgs0,'RMSE'),3);
        npglwRMSE = mean(Cell.getField(  npglw,'RMSE'),3);
        npgslwRMSE = mean(Cell.getField( npgslw,'RMSE'),3);

        npgTime = mean(Cell.getField(    npg,'time'),3);
        npgsTime = mean(Cell.getField(   npgs,'time'),3);
        npg0Time = mean(Cell.getField(   npg0,'time'),3);
        npgs0Time = mean(Cell.getField(  npgs0,'time'),3);
        npglwTime = mean(Cell.getField(  npglw,'time'),3);
        npgslwTime = mean(Cell.getField( npgslw,'time'),3);

        for k=1:4
            for i=1:length(m)
                gnetRMSE(i,k)=min(gnet{i,k}.RMSE(:));
            end
        end
        gnetRMSE=mean(gnetRMSE,2);

        npgsTranRMSE   = zeros(length(m),3);
        npgs0TranRMSE  = zeros(length(m),3);
        npgslwTranRMSE = zeros(length(m),3);
        for k=1:4;
            for i=1:length(npgsRMSE(:,1))
                for j=1:length(npgsRMSE(1,:))
                    npgsTranRMSE(i,j)   = npgsTranRMSE(i,j)   + 1/4*rmseTruncate(  npgs{i,j,k});
                    npgs0TranRMSE(i,j)  = npgs0TranRMSE(i,j)  + 1/4*rmseTruncate( npgs0{i,j,k});
                    npgslwTranRMSE(i,j) = npgslwTranRMSE(i,j) + 1/4*rmseTruncate(npgslw{i,j,k});
                end
            end
        end

        [r,c1]=find(       npgRMSE== repmat(min(       npgRMSE,[],2),1,size(npgRMSE,2))); [r,idx1]=sort(r);
        [r,c2]=find(      npgsRMSE== repmat(min(      npgsRMSE,[],2),1,size(npgRMSE,2))); [r,idx2]=sort(r);
        [r,c3]=find(      npg0RMSE== repmat(min(      npg0RMSE,[],2),1,size(npgRMSE,2))); [r,idx3]=sort(r);
        [r,c4]=find(     npgs0RMSE== repmat(min(     npgs0RMSE,[],2),1,size(npgRMSE,2))); [r,idx4]=sort(r);
        [r,c5]=find(     npglwRMSE== repmat(min(     npglwRMSE,[],2),1,size(npgRMSE,2))); [r,idx5]=sort(r);
        [r,c6]=find(    npgslwRMSE== repmat(min(    npgslwRMSE,[],2),1,size(npgRMSE,2))); [r,idx6]=sort(r);
        [r,c7]=find(  npgsTranRMSE== repmat(min(  npgsTranRMSE,[],2),1,size(npgRMSE,2))); [r,idx7]=sort(r);
        [r,c8]=find( npgs0TranRMSE== repmat(min( npgs0TranRMSE,[],2),1,size(npgRMSE,2))); [r,idx8]=sort(r);
        [r,c9]=find(npgslwTranRMSE== repmat(min(npgslwTranRMSE,[],2),1,size(npgRMSE,2))); [r,idx9]=sort(r);

        disp([c1(idx1), c2(idx2), c3(idx3), c4(idx4)  c5(idx5), c6(idx6) c7(idx7) c8(idx8) c9(idx9)]);

        figure;
        semilogy(m,       npgRMSE(:,2),'r-*'); hold on;
        plot(    m,      npgsRMSE(:,2),'r-p');
        plot(    m,      npg0RMSE(:,2),'g^-');
        plot(    m,     npgs0RMSE(:,2),'gh-');
        plot(    m,     npglwRMSE(:,2),'b<');
        plot(    m,    npgslwRMSE(:,2),'bs-');
        plot(    m,  npgsTranRMSE(:,2),'co-');
        plot(    m, npgs0TranRMSE(:,2),'c.-');
        plot(    m,npgslwTranRMSE(:,2),'ks-');
        legend('npg','npgs','npg0','npgs0','npglw','npgslw','npgsTran','npg0Tran','npglwTran');

        figure;
        semilogy(m,       npgTime(:,2),'r-*'); hold on;
        plot(    m,      npgsTime(:,2),'r-p');
        plot(    m,      npg0Time(:,2),'g^-');
        plot(    m,     npgs0Time(:,2),'gh-');
        plot(    m,     npglwTime(:,2),'b<');
        plot(    m,    npgslwTime(:,2),'bs-');
        legend('npg','npgs','npg0','npgs0','npglw','npgslw');

        figure;
        semilogy(m,mean(       npgRMSE,2),'r-*'); hold on;
        plot(    m,mean(      npgsRMSE,2),'r-p');
        plot(    m,mean(      npg0RMSE,2),'g^-');
        plot(    m,mean(     npgs0RMSE,2),'gh-');
        plot(    m,mean(     npglwRMSE,2),'b<');
        plot(    m,mean(    npgslwRMSE,2),'bs-');
        plot(    m,mean(  npgsTranRMSE,2),'co-');
        plot(    m,mean( npgs0TranRMSE,2),'c.-');
        plot(    m,mean(npgslwTranRMSE,2),'ks-');
        legend('npg','npgs','npg0','npgs0','npglw','npgslw','npgsTran','npg0Tran','npglwTran');

end
end

