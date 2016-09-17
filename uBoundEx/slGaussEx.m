function slGaussEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%
%          Skyline Gaussian Linear example, no background noise
%           Vary the number of measurements, with continuation


if(~exist('op','var')) op='run'; end
switch lower(op)
    case 'run'
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear -regexp '(?i)opt'
        filename = [mfilename '.mat'];

        OPT.maxItr=5e4; OPT.thresh=1e-6; OPT.debugLevel=1;
        m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
        for k=1:1
            for i=1:length(m);
                fprintf('%s, i=%d, k=%d\n','slGaussBound',i,k);

                % the following are under nonnegativity constraints

                OPT.m=m(i); OPT.snr=inf;
                [y,Phi,Phit,Psi,Psit,OPT,~,invEAAt]=loadLinear(OPT,k*100+i);
                initSig = Phit(invEAAt*y);

                p=length(Phit(y));
                PsiM=Utils.getMat(Psi,length(Psit(OPT.trueAlpha)));
                Phity=Phit(y);

                cvx_begin
                    variable a(p)
                    minimize( norm( PsiM'*(Phity+a), inf) )
                    subject to
                        a>=0
                cvx_end

                u_1(i)=cvx_optval;
                
                Pncx=@(x) min(x,0);
                u_2(i)=uBound(Psi,Psit,Pncx,zeros(p,1),-Phit(y));

                opt=OPT;
                initSig=zeros(size(opt.trueAlpha)); ur=u_1(i)*100; ul=0;
                ur_rmse=0; ul_rmse=0;

                while(ur-ul>1e-5)
                    fprintf('%10g <-> %10g\n',ul,ur);
                    fprintf('%10g <-> %10g\n',ul_rmse,ur_rmse);
                    opt.u=(ur+ul)/2; opt.thresh=1e-9;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.PNPG(Phi,Phit,Psi,Psit,y,initSig,opt);
                    rmse=norm(out.alpha)
                    if(rmse<=eps)
                        ur=opt.u; ur_rmse=rmse;
                    else
                        ul=opt.u; ul_rmse=rmse;
                    end
                end
                u_rmse(i)=ur_rmse;
                u_3(i)=ur;

                % the following are under sparsity regularization only
                %
                cvx_begin
                    variable a(p)
                    minimize( norm( PsiM'*(Phity+a), inf) )
                    subject to
                        a>=0
                cvx_end

                u_4(i)=norm( PsiM'*(Phity), inf);
                
                Pncx=@(x) x*0;
                u_5(i)=uBound(Psi,Psit,Pncx,zeros(p,1),-Phit(y));

                opt=OPT;
                initSig=zeros(size(opt.trueAlpha)); ur=u_1(i)*100; ul=0;
                ur_rmse=0; ul_rmse=0;

                while(ur-ul>1e-5)
                    fprintf('%10g <-> %10g\n',ul,ur);
                    fprintf('%10g <-> %10g\n',ul_rmse,ur_rmse);
                    opt.u=(ur+ul)/2; opt.thresh=1e-9;
                    fprintf('u=%g\n',opt.u);
                    out=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt);
                    rmse=norm(out.alpha)
                    if(rmse<=eps)
                        ur=opt.u; ur_rmse=rmse;
                    else
                        ul=opt.u; ul_rmse=rmse;
                    end
                end
                u_rmse(i)=ur_rmse;
                u_6(i)=ur;

                mysave;
            end;
        end

    case 'plot'
        load([mfilename '.mat']);

        m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200

        forSave=[m; u_1; u_2; u_3; u_4; u_5; u_6]';

        figure;
        plot(m,u_1,'b^-'); hold on;
        plot(m,u_2,'gh-'); hold on;
        plot(m,u_3,'bs--');
        plot(m,u_4,'r*-');
        plot(m,u_5,'gp-');
        plot(m,u_6,'ro--');

        rowLabels={'$N$','theoretical','empirical','theoretical','empirical'};
        matrix2latex(forSave(:,[1 2 4 5 7]), 'slBound.tex', 'columnLabels', rowLabels,...
            'alignment', 'r', 'format', '%-6.2f', 'size', 'small');
        save('slBound.data','forSave','-ascii');
end
end

