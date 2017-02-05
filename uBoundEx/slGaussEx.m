function slGaussEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (gurenliang@gmail.com)
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

    OPT.maxItr=1; OPT.debugLevel=0; OPT.thresh=1e-9;
    %m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    % The following corresponding to 30dB 20dB ... -30dB
    snr = [1e3 100 10 1 0.1 0.01 1e-3];
    for k=1:1
      [y,Phi,Phit,Psi,Psit,OPT,~,invEAAt]=loadLinear(OPT,k*100);
      p=length(OPT.trueX);
      PsiM=Utils.getMat(Psi,length(Psit(OPT.trueX)));
      v = randn(OPT.m,1);

      for i=1:length(snr);
        fprintf('%s, i=%d, k=%d\n','slGaussBound',i,k);
        yy = Phi(OPT.trueX)+v*(norm(y)/sqrt(snr(i)*OPT.m));

        Phity=Phit(yy);

        % the following are under sparsity regularization +
        % nonnegativity constraints
%       cvx_begin
%       variable a(p)
%       minimize( norm( PsiM'*(Phity+a), inf) )
%       subject to
%       a>=0
%       cvx_end
%       u_1(i)=cvx_optval;

%       Pncx=@(x) min(x,0);
%       u_2(i)=uBound(Psi,Psit,Pncx,zeros(p,1),-Phity);

%       opt=OPT; opt.maxPossibleInnerItr=1e4;
%       func=@(init,optt) Wrapper.PG(Phi,Phit,Psi,Psit,yy,init,optt);
%       cond=@(x) relativeDif(x,0);
%       u_3(i)=bisection(opt,func,0,u_1(i)*100)

%       % the following are under sparsity regularization only
%       u_4(i)=norm( PsiM'*(Phity), inf);

%       Pncx=@(x) x*0;
%       u_5(i)=uBound(Psi,Psit,Pncx,zeros(p,1),-Phity);

%       opt=OPT;
%       func=@(init,optt) Wrapper.NPGs(Phi,Phit,Psi,Psit,yy,init,optt);
%       cond=@(x) relativeDif(x,0);
%       u_6(i)=bisection(opt,func,cond,0,u_4(i)*100);

%       % following is the 1d TV regularization
        x0=ones(p,1)*sum(Phity)/sqrNorm(Phi(ones(p,1)));
        x0=max(x0,0);  % x0 has to be nonnegative
        g=Phit(Phi(x0)-yy);
        u_7(i)=norm(cumsum(g),inf);

        if(x0(1)>0)
          Pncx=@(x) x*0;
        else
          Pncx=@(x) min(x,0);
        end

        %??? set opt with step size

%       u_8(i)=uBound([],[],'l1',Pncx,x0,g);

        opt=OPT;
        C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
        tvType='l1';
        proximal=tvProximal(tvType,C.prox,'pnpg');
        NLL=@(x) Utils.linearModel(x,Phi,Phit,y);
        opt.debugLevel=2; opt.maxItr=1e2; opt.maxInnerItr=1e4;
        opt.trueX=x0; opt.outLevel=2;
        func=@(u) pnpg(NLL,proximal,x0,setfield(opt,'u',u));
        cond=@(x) relativeDif(x,mean(x));
        u_A(i)=bisection(func,cond,0,u_7(i)*1.2,1e-6);

        keyboard


        opt=[]; opt.debugLevel=2; opt.verbose=1e3;
        iso=tvProximal('iso',@(x)max(0,x),[],opt);
        beta=1/OPT.L;
        func=@(u) iso.prox(x0-beta*g,u*beta,1e-11,1e5,[]);
        cond=@(x) relativeDif(x,mean(x));
        u_9(i)=bisection(func,cond,0,u_7(i)*1.2, 10^-6);
 
        mysave;
        keyboard
      end;
    end

  case 'plot'
    load([mfilename '.mat']);

    m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    snr = [1e3 100 10 1 0.1 0.01 1e-3];

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

function x = A(p)
  [I,J]=size(p);
  p(I,:)=0;
  x=[p(1,:); p(2:I,:)-p(1:I-1,:)];
end

function p = At(x)
  [I,J]=size(x);
  p=[x(1:I-1,:)-x(2:I,:);zeros(1,J)];
end

