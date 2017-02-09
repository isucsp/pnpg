function tv_Bound(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upper-Bounding the Regularization Constant for Sparse Signal Reconstruction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (gurenliang@gmail.com)
%
%   PET example, with background noise b
%   Vary the total counts of the measurements, with continuation

if(~exist('op','var')) op='run'; end

switch lower(op)
case 'run'
    filename = [mfilename '.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear -regexp '(?i)opt'
    filename = [mfilename '.mat'];
    C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
    OPT.maxItr=1e4; OPT.thresh=1e-6; OPT.debugLevel=1; OPT.noiseType='poisson';

    K=1;
    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    for k=1:K
        for i=1:length(count)
            fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,1,k);
            OPT.mask  =[];
            [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
            NLL=@(x) Utils.poissonModel(x,Phi,Phit,y,OPT.bb);

            [x0s,g]=Utils.poissonModelConstEst(Phi,Phit,y,OPT.bb,1e-16);
            if(x0s>0) Pncx=@(x) x*0; else Pncx=@(x) min(x,0); end

%           tvType='l1';
%           u_1(i)=uBound([],[],tvType,Pncx,x0s*ones(size(g)),g);
%           mysave;

%           initSig=ones(size(OPT.trueX))*x0s;
%           PROXOPT=[]; PROXOPT.debugLevel=0; PROXOPT.verbose=1e3;
%           proximal=tvProximal(tvType,C.prox,[],PROXOPT);
%           cond=@(x) relativeDif(x,mean(x(:)));
%           %beta=1/OPT.L;
%           %func=@(u) proximal.prox(x0-beta*g,u*beta,1e-11,1e5,[]);
%           opt=OPT; opt.debugLevel=2; opt.maxInnerItr=1e4;
%           opt.errorType=-1; opt.computError=@(x) relativeDif(x,mean(x(:)));
%           func=@(u) pnpg(NLL,proximal,x0s*ones(size(g)),setfield(opt,'u',u));
%           u_2(i)=bisection(func,cond,u_1(i)/2,u_1(i)*2, 10^-6);

%           tvType='iso';
%           u_4(i)=uBound([],[],tvType,Pncx,x0s*ones(size(g)),g);
%           mysave;

%           initSig=ones(size(OPT.trueX))*x0s;
%           PROXOPT=[]; PROXOPT.debugLevel=0; PROXOPT.verbose=1e3;
%           proximal=tvProximal(tvType,C.prox,[],PROXOPT);
%           cond=@(x) relativeDif(x,mean(x(:)));
%           %beta=1/OPT.L;
%           %func=@(u) proximal.prox(x0-beta*g,u*beta,1e-11,1e5,[]);
%           opt=OPT; opt.debugLevel=2; opt.maxInnerItr=1e4;
%           opt.errorType=-1; opt.computError=@(x) relativeDif(x,mean(x(:)));
%           func=@(u) pnpg(NLL,proximal,x0s*ones(size(g)),setfield(opt,'u',u));
%           u_5(i)=bisection(func,cond,u_4(i)/2,u_4(i)*2, 10^-6);

%           mysave

            Pncx=@(x) min(x,0);
            u_7(i)=uBound(Psi,Psit,'wav',Pncx,zeros(size(OPT.trueX)),...
                Phit(1-y./OPT.bb));
             
            PROXOPT.Lip=@(u)u^2; PROXOPT.initStep='fixed';
            PROXOPT.adaptiveStep=false; PROXOPT.backtracking=false;
            PROXOPT.debugLevel=2; PROXOPT.verbose=1e3;
            proximal=sparseProximal(Psi,Psit,C.prox,'pnpg',PROXOPT);
            cond=@(x) norm(x)/length(x);
            %beta=1/OPT.L;
            beta=1;  % a step size is not needed sine x*=0
            func=@(u) proximal.prox(-beta*Phit(1-y./OPT.bb), u*beta,1e-9,1e5,[]);
            %opt=OPT; opt.maxPossibleInnerItr=1e4; opt.trueX=opt.trueX*0;
            %func=@(u) pg(NLL,proximal,x0*0,setfield(opt,'u',u));
            u_8(i)=bisection(func,cond,u_7(i)/2,u_7(i)*2,1e-6);

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

