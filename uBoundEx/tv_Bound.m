function tv_Bound(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%
%
%                  PET example, with background noise b
%      Vary the total counts of the measurements, with continuation

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear -regexp '(?i)opt'
        filename = [mfilename '.mat'];
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1; opt.noiseType='poisson';
        opt.minItr=30;
        opt.mask  =[];

        K=1;
        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        for k=1:K
            for i=4:length(count)
                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,1,k);
                [y,Phi,Phit,Psi,Psit,fbpfunc,opt]=loadPET(count(i),opt);

                [x0s,g]=Utils.poissonModelConstEst(Phi,Phit,y,opt.bb,1e-16);
                g=reshape(g,sqrt(length(g(:))),[]);

                %u_maxANI(i,k)=TV.upperBoundU(maskFunc(g,opt.mask));
                %u_maxISO(i,k)=sqrt(2)*u_maxANI(i,k);
                %u_maxANI_dual(i,k)=TV.upperBoundU_dual(maskFunc(g,opt.mask));
                %u_maxISO_dual(i,k)=sqrt(2)*u_maxANI_dual(i,k);
                %u_maxANI_admm(i,k)=TV.upperBoundU_admm(maskFunc(g,opt.mask));
                %u_maxANI_adm2(i,k)=TV.upperBoundU_admm2(maskFunc(g,opt.mask));
                %keyboard
                u_admm(i,k)=TV.upperBoundU_admm3(g,x0s*ones(size(g)));

                continue;
                initSig=ones(size(opt.trueAlpha))*x0s;

                rmse=0; opt.u=u_maxANI(i,k); opt.proximal='tvl1';
                while(rmse==0)
                    opt.u = 0.9*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                opt.u=opt.u/0.9; rmse=0;
                while(rmse==0)
                    opt.u = 0.99*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                u_trueANI(i,k)=opt.u/0.99;

                rmse=0; opt.u=u_maxISO(i,k); opt.proximal='tviso';
                while(rmse==0)
                    opt.u = 0.9*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                opt.u=opt.u/0.9; rmse=0;
                while(rmse==0)
                    opt.u = 0.99*opt.u;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPG(Phi,Phit,[],[],y,initSig,opt);
                    rmse=norm(out.alpha-initSig);
                end
                u_trueISO(i,k)=opt.u/0.99;
            end
        end

        figure;
        loglog(count,u_maxANI,'b^-'); hold on;
        loglog(count,u_maxANI_dual,'gh-'); hold on;
        loglog(count,u_trueANI,'bs--');
        loglog(count,u_maxISO,'r*-');
        loglog(count,u_maxISO_dual,'gp-');
        loglog(count,u_trueISO,'ro--');
        h=legend('U_0','empirical anisotropic U','sqrt(2)U_0','empirical isotropic U');
        set(h,'interpreter','latex');

        forSave=[count(:) u_maxANI, u_trueANI, u_maxISO, u_trueISO];
        save('bound4U.data','forSave','-ascii');

        save(filename);
end

end

