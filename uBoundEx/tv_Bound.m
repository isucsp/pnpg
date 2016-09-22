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
        OPT.maxItr=1e4; OPT.thresh=1e-6; OPT.debugLevel=1; OPT.noiseType='poisson';

        K=1;
        count = [1e4 1e5 1e6 1e7 1e8 1e9];
        for k=1:K
            for i=1:length(count)
                fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,1,k);
        %       OPT.mask  =[];
        %       [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);

        %       [x0s,g]=Utils.poissonModelConstEst(Phi,Phit,y,OPT.bb,1e-16);
        %       g=reshape(g,sqrt(length(g(:))),[]);

        %       u_1(i)=TV.upperBoundU_admm3(g,x0s*ones(size(g)));

        %       initSig=ones(size(OPT.trueAlpha))*x0s;
        %       ur=u_1(i)*10; ul=ur/100; ur_rmse=0; ul_rmse=0; opt=OPT;
        %       opt.proximal='tvl1'; opt.maxItr=13;
        %       while(ur-ul>1e-5*ur)
        %           fprintf('%10g <-> %10g\n',ul,ur);
        %           fprintf('%10g <-> %10g\n',ul_rmse,ur_rmse);
        %           opt.u=(ur+ul)/2; opt.thresh=1e-9;
        %           fprintf('u=%g\n',opt.u);
        %           out=Wrapper.PNPG(Phi,Phit,[],[],y,initSig,opt);
        %           rmse=norm(out.alpha-initSig)
        %           if(rmse<=eps)
        %               ur=opt.u; ur_rmse=rmse;
        %           else
        %               ul=opt.u; ul_rmse=rmse;
        %           end
        %       end
        %       u_2(i)=ur;
        %       u_2rmse(i)=ur_rmse;

                if(isfield(OPT,'mask')) OPT=rmfield(OPT,'mask'); end;
                [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
                Pncx=@(x) min(x,0);
                %u_3(i)=uBound(Psi,Psit,Pncx,OPT.trueAlpha*0,Phit(1-y./OPT.bb(:)));
                u_3(i)=1e3;

                ur=u_3(i)*10; ul=ur/100; ur_rmse=0; ul_rmse=0; opt=OPT;
                initSig=opt.trueAlpha*0; opt.maxItr=13;
                opt.proximal='wvltADMM';
                while(ur-ul>1e-5*ur)
                    fprintf('%10g <-> %10g\n',ul,ur);
                    fprintf('%10g <-> %10g\n',ul_rmse,ur_rmse);
                    opt.u=(ur+ul)/2; opt.thresh=1e-9;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.PNPG(Phi,Phit,Psi,Psit,y,initSig,opt);
                    rmse=norm(out.alpha-initSig)
                    if(rmse<=eps)
                        ur=opt.u; ur_rmse=rmse;
                    else
                        ul=opt.u; ul_rmse=rmse;
                    end
                end
                u_4(i)=ur;
                u_4rmse(i)=ur_rmse;

                if(i==1) keyboard; end
                mysave;
            end
        end

    case 'plot'
        load([mfilename '.mat']);

        count = [1e4 1e5 1e6 1e7 1e8 1e9];

        figure;
        loglog(count,u_1,'b^-'); hold on;
        loglog(count,u_2,'bs--');
        loglog(count,u_1*sqrt(2),'gh-'); hold on;
        loglog(count,u_3,'r*-');
        h=legend('$U_0$','empirical anisotropic $U$','$\sqrt{2}U_0$','empirical isotropic $U$');
        set(h,'interpreter','latex');

        forSave=[count(:), u_1', u_2', u_3'];
        save('bound4U.data','forSave','-ascii');
        save(filename);
end
end
