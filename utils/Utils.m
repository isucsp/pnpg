classdef Utils < handle
methods(Static)
function out = getNthOutPut(f,x,n)
    switch(n)
        case 2
            [~,out]=f(x);
        case 3
            [~,~,out]=f(x);
        case 4
            [~,~,~,out]=f(x);
    end
end
function y = softThresh(x,thresh)
    thresh=abs(thresh);
    y=zeros(size(x));
    idx=(x>thresh);
    y(idx)=x(idx)-thresh;
    idx=(x<-thresh);
    y(idx)=x(idx)+thresh;
end
function mat = getMat(func,ncol)
    mat = zeros(length(func(zeros(ncol,1))),ncol);
    x = zeros(ncol,1);
    str='';
    fprintf('Function Handle to Matrix Recover: ');
    for j=1:ncol
        strlen=length(str); str=sprintf('%d-th column ...',j); fprintf([repmat('\b',1,strlen) '%s'],str);
        x(j)=1; mat(:,j)=func(x); x(j)=0;
    end
    fprintf('  done! \n');
end
function [f,g,h] = augLag(x,z)
    g=x-z;
    f=sqrNorm(x-z)/2;
    if(nargout>=3)
        h = ones(size(x));
    end
end
function varargout = linearModel(alpha,Phi,Phit,y,b)
    PhiAlpha=Phi(alpha);
    if(exist('b','var')) PhiAlpha=PhiAlpha+b; end
    switch(nargout)
        case 2
           [varargout{1},g]=Utils.augLag(PhiAlpha,y);
           varargout{2}=Phit(g);
        case 3
           [varargout{1},g]=Utils.augLag(PhiAlpha,y);
           varargout{2}= Phit(g);
           varargout{3}= @(x,opt) hessian(x,opt);
       otherwise
           varargout{1}=Utils.augLag(PhiAlpha,y);
    end
    function hh = hessian(x,opt,h)
        if(opt==1)
            hh = Phit(Phi(x));
        else
            hh = sqrNorm(Phi(x));
        end
    end
end
function [f,g,h] = nonnegLogBarrier(alpha)
    %if(any(alpha(:)<=0)) f=eps^-1; alpha(alpha<=0)=eps;
    %else f=-sum(log(alpha(:)));
    %end
    %if(nargout>1) g = -1./alpha; h=1./alpha.^2; end
    f=-sum(log(alpha(:)));
    if(nargout>1)
        g=-1./alpha;
        h=1./(alpha.^2);
    end
end
function [f,g,h] = nonnegPen(alpha)
    temp=(alpha<0);
    f=alpha(temp)'*alpha(temp);
    if(nargout>=2)
        g=zeros(size(alpha));
        g(temp)=2*alpha(temp);
        if(nargout>=3)
            h = @(x,opt) hessian(x,opt);
        end
    end
    function hh = hessian(x,opt)
        if(opt==1)
            hh = zeros(size(x));
            hh(temp,:) = x(temp,:)*2;
        else
            x = reshape(x,length(temp(:)),[]);
            y = x(temp,:);
            hh = y'*y*2;
        end
    end
end

function [f,g,h]=barrierIe(Ie)
    %if(any(Ie)<=0)
    %    Ie(Ie<=0)=eps; f=eps^-1;
    %    if(1-sum(Ie)<=0) Ie=Ie*(1-eps)/sum(Ie); end
    %else
    %    if(1-sum(Ie)<=0) Ie=Ie*(1-eps)/sum(Ie); f=eps^-1;
    %    else f=-sum(log(Ie))-log(1-sum(Ie)); end
    %end
    f=-sum(log(Ie))-log(1-sum(Ie));
    if(nargout>1)
        g=1/(1-sum(Ie))-1./Ie;
        h=1/(1-sum(Ie))^2+diag(1./(Ie.^2));
    end
end

function [f,g,h] = lustigL1(alpha,xi,Psi,Psit)
    s=Psit(alpha);
    sqrtSSqrMu=sqrt(s.^2+xi);
    f=sum(sqrtSSqrMu);
    if(nargout>=2)
        g=Psi(s./sqrtSSqrMu);
        if(nargout>=3)
            weight = xi./(sqrtSSqrMu.^3);
            h = @(x,opt) hessian(weight,x,opt);
        end
    end
    function hh = hessian(weight,x,opt)
        y = Psit(x);
        if(opt==1)
            hh = Psi(weight.*y);
        else
            hh = y'*(weight.*y);
        end
    end
end

function [f,g,h] = huber(alpha,mu,Psi,Psit)
    s=Psit(alpha);
    idx = abs(s)<mu;
    temp = abs(s)-mu/2;
    temp(idx) = s(idx).^2/2/mu;
    f = sum(temp);
    if(nargout>=2)
        temp = ones(size(s));
        temp(idx) = s(idx)/mu;
        g=Psi(temp);
        if(nargout>=3)
            h = @(x,opt) hessian(x,opt);
        end
    end
    function hh = hessian(x,opt)
        y = Psit(x);
        if(opt==1)
            y(idx) = y(idx)/mu; y(~idx) = 0; hh=Psi(y);
        else
            y = y(idx); hh = y'*y/mu;
        end
    end
end

% The Poisson Generalized Linear Model with identity log link: Ey=Φα
function [f,g,h] = poissonModelAppr(alpha,Phi,Phit,y,b,EPS)
    if(~exist('b','var')) b=0; end
    if(exist('EPS','var')) eps=EPS; else eps=1e-10; end

    PhiAlpha=Phi(alpha)+b;
    nzpx=(PhiAlpha>eps);
    nzy=(y>0);
    f=sum(PhiAlpha(:)-y);
    f=f-innerProd(y(nzpx & nzy),log(PhiAlpha(nzpx & nzy)./y(nzpx & nzy)));
    f=f-innerProd(y(~nzpx & nzy),...
        log(eps)-1.5+2*PhiAlpha(~nzpx & nzy)/eps-(PhiAlpha(~nzpx & nzy).^2)/4/eps^2-log(y(~nzpx & nzy)));

    if(nargout>=2)
        t=PhiAlpha;
        t(~nzpx)=eps./(2-t(~nzpx)/eps);
        t=1-y./t;
        g=Phit( t );
        if(nargout>=3)
            weight=max(PhiAlpha,eps);
            weight=y./(weight.^2);
            h=@(x,opt) hessian(weight,x,opt);
        end
    end
    function hh=hessian(weight,x,opt)
        z = Phi(x);
        if(opt==1)
            hh = Phit(weight.*z);
        else
            hh = innerProd(z,weight.*z);
        end
    end
end

function [f,g,h] = poissonModel(alpha,Phi,Phit,y,b,EPS)
    % modeling the Poisson measurement with linear model
    %     y ~ Poisson ( Φα + b )
    %  f(α)=1'*(Φα+b) - y'*ln(Φα+b)
    %  g(α)=Φ'*[ 1-y./(Φα+b) ]
    %  default eps is zero, suggest to just use default value
    if(exist('EPS','var')) eps=EPS; else eps=0; end
    if(~exist('b','var')) b=0; end

    PhiAlpha=Phi(alpha)+b;
    nzy=(y~=0);
    nzpx=(PhiAlpha~=0);
    f=sum(PhiAlpha(:)-y(:))-innerProd(y(nzy),log(PhiAlpha(nzy)./y(nzy)));

    if(nargout>=2)
        t=ones(size(y));
        t(nzy)=1-y(nzy)./(PhiAlpha(nzy)+eps);
        g=Phit( t );
        if(nargout>=3)
            weight=zeros(size(y));
            weight(nzy)=y(nzy)./((PhiAlpha(nzy)+eps).^2);
            h=@(x,opt) hessian(weight,x,opt);
        end
    end
    function hh=hessian(weight,x,opt)
        z = Phi(x);
        if(opt==1)
            hh = Phit(weight.*z);
        else
            hh = innerProd(z,weight.*z);
        end
    end
end
function [f,g,h] = poisson(x,y,EPS)
    % modeling the Poisson measurement with linear model
    %     y ~ Poisson ( x )
    %  f(α)=1'*(x-y) - y'*ln(x./y)
    %  g(α)=[ 1-y./x ]
    %  h(α)=[ y./x.^2 ]

    if(nargout>=1)
        nzy=(y~=0);  %nzx=(x~=0);
        f=sum(x-y)-innerProd(y(nzy),log(x(nzy)./y(nzy)));
        if(nargout>=2)
            g=ones(size(y));
            g(nzy)=1-y(nzy)./x(nzy);
            if(nargout>=3)
                h=zeros(size(y));
                h(nzy)=y(nzy)./(x(nzy).^2);
            end
        end
    end
end
function varargout = glm(alpha,y,model,Phi,Phit,b)
    PhiAlpha=Phi(alpha);
    if(exist('b','var')) PhiAlpha=PhiAlpha+b; end
    switch(nargout)
        case 1
           varargout{1}=model(PhiAlpha,y);
        case 2
           [varargout{1},g]=model(PhiAlpha,y);
           varargout{2}=Phit(g);
        case 3
           [varargout{1},g,h]=model(PhiAlpha,y);
           varargout{2}= Phit(g);
           varargout{3}= @(x,opt) hessian(x,opt,h);
    end
    function hh = hessian(x,opt,h)
        yy = Phi(x);
        if(opt==1)
            hh = Phit(h.*yy);
        else
            hh = yy'*(h.*yy);
        end
    end
end
function [a,b] = poissonModelConstEst(Phi,Phit,y,b,epsilon)
    % Estimate the MLE given the signal is a constant
    if(~exist('b','var') || isempty(b)) b=0; end

    one = ones(size(Phit(y)));
    PhiOne=Phi(one);
    OnePhiOne=sum(reshape(PhiOne,[],1));
    yPhiOne=y.*PhiOne;
    g=@(x) OnePhiOne-sum(reshape(yPhiOne./(PhiOne*x+b),[],1));

    high=1;
    if(g(high)>0)
        low=high/10;
        while(g(low)>0)
            high=low;
            low=low/10;
        end
    elseif(g(high)<0)
        while(g(high)<0)
            low=high;
            high=high*10;
        end
    else
        low=high;
    end
    while(high-low>epsilon)
        mid=(high+low)/2;
        if(g(mid)>0)
            high=mid;
        elseif(g(mid)<0)
            low=mid;
        else
            high=mid; low=mid;
        end
    end
    a=high;
    [~,b]=Utils.poissonModel(one*a,Phi,Phit,y,b);
end

% The following model assumes the y=I_0 * exp(-Φα), where I_0 is known
function [f,g,h] = poissonModelLogLink0(alpha,Phi,Phit,y,I0)
    PhiAlpha=Phi(alpha); weight = exp(log(I0)-PhiAlpha);
    nzy=(y~=0);
    f=sum(weight-y)+innerProd(y(nzy),PhiAlpha(nzy)+log(y(nzy)/I0));
    if(nargout>=2)
        g=Phit(  y(:) - weight  );
        if(nargout>=3)
            h=@(x,opt) hessian(weight,x,opt);
        end
    end
    function hh=hessian(weight,x,opt)
        z = Phi(x);
        if(opt==1)
            hh = Phit(weight.*z);
        else
            hh = innerProd(z,weight.*z);
        end
    end
end

% The following model assumes the y=I_0 * exp(-Φα), where I_0 is unknown and estimated as a function of α
function [f,g,h] = poissonModelLogLink(alpha,Phi,Phit,y)
    PhiAlpha=Phi(alpha); weight = exp(-PhiAlpha);
    sumy=sum(y); sumWeight=sum(weight);
    nzy=(y~=0);
    f=innerProd(y(nzy), PhiAlpha(nzy)+log(y(nzy)*sumWeight/sumy));
    % f=sumy*log(sumWeight)+innerProd(y,PhiAlpha);
    if(nargout>=2)
        g=Phit(  y(:) - weight*sumy/sumWeight  );
        if(nargout>=3)
            h=@(x,opt) hessian(x,opt);
        end
    end
    function hh = hessian(x,opt)
        z = Phi(x);
        if(opt==1)
            hh = Phit(weight.*z-weight*innerProd(weight,z)/sumWeight)*sumy/sumWeight;
        else
            w = innerProd(weight,z);
            hh = (innerProd(z,weight.*z)-w*w/sumWeight)*sumy/sumWeight;
        end
    end
end

% logistic model where b is known
function [f,g,h] = logisticModel(x,Phi,Phit,y,b)
    if(nargin<5) b=zeros(size(y)); end
    PhiX=Phi(x)+b;
    expPhiX=exp(PhiX);
    f=mean(-y.*PhiX+log(1+expPhiX));
    if(nargout>=2)
        t=-y+expPhiX./(1+expPhiX);
        g=Phit( t )/length(y);
        if(nargout>=3)
            weight = expPhiX./((1+expPhiX).^2)/length(y);
            h=@(x,opt) hessian(weight,x,opt);
        end
    end
    function hh=hessian(weight,x,opt)
        z = Phi(x);
        if(opt==1)
            hh = Phit(weight.*z);
        else
            hh = innerProd(z,weight.*z);
        end
    end
end

function [dif,l1norm]=testSparsity(x)
    k=0;
    minRSE=inf; maxRSE=0;
    n=max(size(x));
    if(n==length(x(:)))
        L=floor(log2(n));
    else
        L=floor(log2(min(size(x))));
    end
    for perce=0.7:0.01:0.98
        k=k+1;
        for i=1:min(L,9)
            for j=1:4
                level = i; wave = 2*j;
                [dif(i,j,k),~,l1norm(i,j)]=Utils.getError(x,perce,level,wave);
            end
        end
        minDif=dif(:,:,k); minDif=minDif(:);
        if(any(minDif<=0)) keyboard; end
        minRSE = min( minRSE, min(minDif(minDif>0)));
        maxRSE = max( maxRSE, max(minDif(minDif>0)));
        minDif=min(minDif(minDif>0));
        if(~isscalar(minDif)) keyboard; end
        [i,j]=find(dif(:,:,k)==minDif);
        fprintf('level=%d,daub=%d  :  %g%% of signal achieves %g%%\n',i,2*j,(1-perce)*100,dif(i,j,k)*100);
    end
    [i,j]=find(l1norm==min(l1norm(l1norm>0)));
    fprintf('level=%d,daub=%d  :  |x|_1 of signal achieves %g\n',i,2*j,l1norm(i,j));

    k=1;
    for rse = logspace(log10(1e-5), log10(0.5), 20)
        for i=1:min(L,9)
            for j=1:4
                level = i; wave = 2*j;
                low=0; high=1;
                while(high-low>1e-10)
                    perce=(low+high)/2;
                    [dif1(i,j,k),coef]=Utils.getError(x,perce,level,wave);
                    if(dif1(i,j,k)<=rse) low=perce; else high=perce; end
                end
                l1norm1(i,j,k)=sum(abs(coef(1:ceil((1-perce)*length(coef)))));
            end
        end
        minNorm=min(reshape(l1norm1(:,:,k),[],1));
        [i,j]=find(l1norm1(:,:,k)==minNorm);
        fprintf('level=%d,daub=%d  :  |x|_1=%g for signal achieves %g%%\n',i,2*j,minNorm,dif1(i,j,k)*100);

        k=k+1;
    end
end

function [dif,coef,l1norm]=getError(x,perce,level,wave)
    h=daubcqf(wave);
    [wav,~] = mdwt(x,h,level);
    [~,idx]=sort(abs(wav(:)));
    coef=wav(idx);
    coef=coef(end:-1:1);
    %figure(2); plot(wav(idx)); hold on;
    wav(idx(1:floor(perce*length(x(:)))))=0;
    recX=midwt(wav,h,level);
    %figure(1); plot(recX); hold on; plot(x,'r');
    dif=sqrNorm(x-recX)/max(sqrNorm(x),1e-10);
    l1norm=sum(abs(coef(:)));
end

function mask = getCircularMask(n,m)
    if(nargin==1) m=n; end
    mask = zeros(n);
    [x,y]=meshgrid(0:n-1);
    idx = sqrt((x-n/2).^2+(y-n/2).^2)<(m/2-1);
    mask(idx)=1;
end
function [Psi,Psit] = getPsiPsit(daub,dwt_L,maskmt,maskkmt)
    if(~exist('maskmt','var') || isempty(maskmt))
        mask.a=@(xx) xx;
        mask.b=@(xx) xx;
    else
        maskIdx=find(maskmt~=0);
        n=size(maskmt);
        mask.a=@(xx) maskFunc(xx,maskIdx);
        mask.b=@(xx) maskFunc(xx,maskIdx,n);
    end
    if(~exist('maskkmt','var') || isempty(maskkmt))
        maskk.a=@(xx) xx;
        maskk.b=@(xx) xx;
    else
        maskkIdx=find(maskkmt~=0);
        m=size(maskkmt);
        maskk.a=@(xx) maskFunc(xx,maskkIdx);
        maskk.b=@(xx) maskFunc(xx,maskkIdx,m);
    end
    wav=daubcqf(daub);

    %Sampling operator
    Psi =@(z) mask.a(midwt(maskk.b(z),wav,dwt_L));
    Psit=@(z) maskk.a(mdwt(mask.b(z),wav,dwt_L));
end
end
end

