function slGaussEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upper-Bounding the Regularization Constant for Sparse Signal Reconstruction %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (gurenliang@gmail.com)
%
%   Skyline Gaussian Linear example, no background noise
%   Vary SNR of measurements, with continuation


if(~exist('op','var')) op='run'; end
switch lower(op)
  case 'run'
    filename = [mfilename '.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear -regexp '(?i)opt'
    filename = [mfilename '.mat'];

    C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
    OPT.maxItr=1; OPT.debugLevel=0; OPT.thresh=1e-9;
    %m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    % The following corresponding to 30dB 20dB ... -30dB
    snr = [ 1e3 100 10 1 0.1 0.01 1e-3 1e-3 1e5 1e4 ];
    for k=1
      [y,Phi,Phit,Psi,Psit,OPT,~,invEAAt]=loadLinear(OPT,k*100);
      p=length(OPT.trueX);
      PsiM=Utils.getMat(Psi,length(Psit(OPT.trueX)));
      v = randn(OPT.m,1);

      for i=8:length(snr);
        fprintf('%s, i=%d, k=%d\n','slGaussBound',i,k);
        yy = Phi(OPT.trueX)+v*(norm(y)/sqrt(snr(i)*OPT.m));
        if(i==8)
            % to get a case when x0<0
            for e=1:100
                v = randn(OPT.m,1);
                yy = Phi(OPT.trueX)+v*(norm(y)/sqrt(snr(i)*OPT.m));
                Phity=Phit(yy);
                x0=ones(p,1)*sum(Phity)/sqrNorm(Phi(ones(p,1)));
                if(x0(1)<0) break; end
            end
        end
        NLL=@(x) Utils.linearModel(x,Phi,Phit,yy);
        Phity=Phit(yy);

        % the following are under sparsity regularization +
        % nonnegativity constraints
        cvx_begin
        variable a(p)
        minimize( norm( PsiM'*(Phity+a), inf) )
        subject to
        a>=0
        cvx_end
        u_1(k,i)=cvx_optval;
        fprintf('u_1=%20.10g\n',u_1(k,i));

        %cvx_begin
        %variable v
        %variable w(size(Phity))
        %variable t(size(Phity))
        %minimize( -v )
        %subject to
        %    -v*Phity+PsiM*w+t==0
        %    norm(w,inf)<=1
        %    t<=0
        %cvx_end
        %u_2(i,k)=-1/cvx_optval;
        Pncx=@(x) min(x,0);
        u_2(k,i)=uBound(Psi,Psit,'wav',Pncx,zeros(p,1),-Phity);
        fprintf('u_2=%20.10g\n',u_2(k,i));

        PROXOPT.Lip=@(u)u^2; PROXOPT.initStep='fixed';
        PROXOPT.adaptiveStep=false; PROXOPT.backtracking=false;
        PROXOPT.debugLevel=2; PROXOPT.verbose=1e3;
        proximal=sparseProximal(Psi,Psit,C.prox,'pnpg',PROXOPT);
        cond=@(x) norm(x)/length(x);
        %beta=1/OPT.L;
        beta=1;  % a step size is not needed sine x*=0
        func=@(u) proximal.prox(beta*Phity,u*beta,1e-10,1e5,[]);
        %opt=OPT; opt.maxPossibleInnerItr=1e4; opt.trueX=opt.trueX*0;
        %func=@(u) pg(NLL,proximal,opt.trueX*0,setfield(opt,'u',u));
        u_3(k,i)=bisection(func,cond,u_1(k,i)/2,u_1(k,i)*2,1e-6);
        fprintf('u_3=%20.10g\n',u_3(k,i));

        % the following are under sparsity regularization only
        u_4(k,i)=norm( PsiM'*(Phity), inf);
        fprintf('u_4=%20.10g\n',u_4(k,i));

        Pncx=@(x) x*0;
        u_5(k,i)=uBound(Psi,Psit,'wav',Pncx,zeros(p,1),-Phity);
        fprintf('u_5=%20.10g\n',u_5(k,i));

        proximal.exact=true;
        proximal.val=@(x) norm(Psit(x),1);
        proximal.prox=@(x,u) Psi(Utils.softThresh(Psit(x),u));
        cond=@(x) norm(x)/length(x);
        func=@(u) Psi(Utils.softThresh(Psit(Phity),u));
        %opt=OPT; opt.maxPossibleInnerItr=1e4; opt.trueX=opt.trueX*0;
        %func=@(u) pg(NLL,proximal,opt.trueX*0,setfield(opt,'u',u));
        u_6(k,i)=bisection(func,cond,u_4(k,i)/2,u_4(k,i)*2,1e-6);
        fprintf('u_6=%20.10g\n',u_6(k,i));
        mysave;

        % following is the 1d TV regularization with R_+ constrains
        x0=ones(p,1)*sum(Phity)/sqrNorm(Phi(ones(p,1)));
        x00(k,i)=x0(1);
        x0=max(x0,0);  % x0 has to be nonnegative
        g=Phit(Phi(x0)-yy);

        if(x0(1)==0)
            cvx_begin
            variable a(p)
            minimize( norm( cumsum(g+a), inf) )
            subject to
                a<=0
            cvx_end
            u_7(k,i)=cvx_optval;
        else
            u_7(k,i)=norm(cumsum(g),inf);
        end
        fprintf('u_7=%20.10g\n',u_7(k,i));

        if(x0(1)>0)
          Pncx=@(x) x*0;
        else
          Pncx=@(x) min(x,0);
        end
        optu_8.maxItr=5e4;
        u_8(k,i)=uBound([],[],'l1',Pncx,x0,g,optu_8);
        fprintf('u_8=%20.10g\n',u_8(k,i));

        PROXOPT=[]; PROXOPT.debugLevel=2; PROXOPT.verbose=1e3;
        proximal=tvProximal('iso',C.prox,[],PROXOPT);
        cond=@(x) relativeDif(x,mean(x));
        beta=1/OPT.L;
        func=@(u) proximal.prox(x0-beta*g,u*beta,1e-11,1e5,[]);
        %opt=OPT; opt.debugLevel=2; opt.maxInnerItr=1e4;
        %opt.errorType=-1; opt.computError=@(x) relativeDif(x,mean(x));
        %func=@(u) pnpg(NLL,proximal,x0,setfield(opt,'u',u));
        if(i==8) thresh=1e-5; else thresh=1e-6; end
        u_9(k,i)=bisection(func,cond,u_7(k,i)/2,u_7(k,i)*1.2, thresh);
        fprintf('u_9=%20.10g\n',u_9(k,i));
        mysave;
         
        % following is the 1d TV regularization without nonnegativity
        % constraints
        x0=ones(p,1)*sum(Phity)/sqrNorm(Phi(ones(p,1)));
        if(x0(1)>0) continue; end
        g=Phit(Phi(x0)-yy);
        u_a(k,i)=norm(cumsum(g),inf);
        fprintf('u_a=%20.10g\n',u_a(k,i));

        Pncx=@(x) x*0;
        optu_b.maxItr=5e4;
        u_b(k,i)=uBound([],[],'l1',Pncx,x0,g,optu_b);
        fprintf('u_b=%20.10g\n',u_b(k,i));

        PROXOPT=[]; PROXOPT.debugLevel=2; PROXOPT.verbose=1e3;
        proximal=tvProximal('iso',@(x)x,[],PROXOPT);
        cond=@(x) relativeDif(x,mean(x));
        beta=1/OPT.L;
        func=@(u) proximal.prox(x0-beta*g,u*beta,1e-11,1e5,[]);
        %opt=OPT; opt.debugLevel=2; opt.maxInnerItr=1e4;
        %opt.errorType=-1; opt.computError=@(x) relativeDif(x,mean(x));
        %func=@(u) pnpg(NLL,proximal,x0,setfield(opt,'u',u));
        if(i==8) thresh=1e-5; else thresh=1e-6; end
        u_c(k,i)=bisection(func,cond,u_a(k,i)/2,u_a(k,i)*1.2, thresh);
        fprintf('u_c=%20.10g\n',u_c(k,i));
        mysave;
      end;
    end

  case 'plot'
    load([mfilename '.mat']);

    m = [ 200, 250, 300, 350, 400, 500, 600, 700, 800]; % should go from 200
    snr = [1e3 100 10 1 0.1 0.01 1e-3 1e-3];

    forSave=[u_1; u_2; u_3; u_4; u_5; u_6; u_7; u_8; u_9; u_a; u_b; u_c; snr]';

    %figure;
    %loglog(snr,u_1,'b^-'); hold on;
    %plot(snr,u_2,'gh-'); hold on;
    %plot(snr,u_3,'bs--');
    %plot(snr,u_4,'r*-');
    %plot(snr,u_5,'gp-');
    %plot(snr,u_6,'ro--');
    %plot(snr,u_7,'r*-');
    %plot(snr,u_8,'gp-');
    %plot(snr,u_9,'ro--');

    rowLabels={'$N$','theoretical','empirical',...
        'theoretical','empirical','theoretical','empirical','theoretical','empirical'};

    forSave(:,13)=10*log10(forSave(:,13));
    matrix2latex(forSave(:,[13 1 3 4 6 7 9 10 12]), 'slBound.tex', 'columnLabels', rowLabels,...
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

