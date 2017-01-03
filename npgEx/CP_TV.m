function out=CP_TV(Phi,Phit,y,approach,tvType,G,xInit,opt)
 
    switch lower(opt.noiseType)
        case 'poisson'
            F_=poissonProximal(y,opt.bb);
        case 'gaussian'
            F_=gaussianProximal(y);
    end

    if(approach==1)
        F=F_;
        K.forward=Phi;
        K.backward=Phit;
        proximal=tvProximal(tvType,G.prox);
        G.iterative=proximal.iterative;
        G.val=@(x)opt.u*proximal.val(x);
        G.prox=@(x,u,varargin) proximal.prox(x,u*opt.u,varargin{:});
    else
        F.iterative=false;
        switch(lower(tvType))
            case 'iso'
                F.val=@(z) F_.val(z{1})+opt.u*isoNorm(z{2},z{3});
                F.proxConj=@(a,u) isoPrjCell(F_.proxConj(a{1},u),a{2},a{3},opt.u);
            case 'l1'
                F.val=@(z) F_.val(z{1})+...
                    opt.u*norm([reshape(z{2},[],1); reshape(z{3},[],1)],1);
                F.proxConj=@(a,u) {F_.proxConj(a{1},u),...
                    min(max(a{2},-opt.u),opt.u),min(max(a{3},-opt.u),opt.u)};
        end

        K.forward=@(x) {Phi(x),TV.Psi_vt(x),TV.Psi_ht(x)};
        K.backward=@(y) Phit(y{1})+TV.Psi_v(y{2})+TV.Psi_h(y{3});
    end
    out = chambollePock(F,G,K,xInit,opt);
    out.F=F;
    out.F_=F_;
    out.K=K;
    out.isoNorm=@isoNorm;

end %function out=CP_TV(Phi,Phit,y,approach,tvType,xInit,opt)

function F = poissonProximal(y,b)
    % modeling the Poisson measurement
    %     y ~ Poisson ( x + b )
    % f(x)=1'*(x+b-y) - y'*ln((x+b)./y)
    % F.val=@f;
    % F.prox solves 0.5*||x-a||_2^2+u*f(x);

    if(any(y<0))
        error(['measurements contain negative components,', ...
            ' thus not qualified as Poisson!']);
    end
    nzy=(y~=0);  %nzx=(x~=0);
    sumY=sum(y(:));
    if(isscalar(b))
        sumB=length(y(:))*b;
        F.val=@(x) sum(x(:))+sumB-sumY-y(nzy).'*log((x(nzy)+b)./y(nzy));
    else
        sumB=sum(b(:));
        F.val=@(x) sum(x(:))+sumB-sumY-y(nzy).'*log((x(nzy)+b(nzy))./y(nzy));
    end
    F.iterative=false;
    F.prox=@(a,u) max(0,0.5*( (a-b-u)+sqrt( (a+b-u).^2+4*u*y ) ));
    %F.proxConj=@(a,u) 0.5*( (a+b*u+1)-sqrt( (a-b*u-1).^2+4*u*(y-b)+4*a.*b ) );
    F.proxConj=@(a,u) a-F.prox(a/u,1/u)*u;

end %function F = poissonProximal(y,b0)

function F = gaussianProximal(y)
    % modeling the Gaussian measurement
    %     y ~ Normal ( x )
    % f(x)=0.5*||x-y||_2^2
    % F.val=@f;
    % F.prox solves 0.5*||x-a||_2^2+u*f(x);

    F.iterative=false;
    F.val=@(x) sum(reshape((x-y).^2,[],1))/2;
    F.prox=@(a,u) (a+u*y)/(1+u);
    F.proxConj=@(a,u) (a-u*y)/(1+u);

end

function o = isoNorm(p,q)

    [I,J]=size(p);
    I=I+1;
    o=sum(reshape(sqrt(p(:,1:J-1).^2+q(1:I-1,:).^2),[],1));
    o=o+norm(p(:,J),1)+norm(q(I,:),1);

end

function pq=isoPrjCell(FopConj,p,q,u)

    [I,J]=size(p);
    I=I+1;
    mag=max(1,sqrt(p(:,1:J-1).^2+q(1:I-1,:).^2)/u);
    p(:,1:J-1)=p(:,1:J-1)./mag; p(:,J)=min(max(p(:,J),-u),u);
    q(1:I-1,:)=q(1:I-1,:)./mag; q(I,:)=min(max(q(I,:),-u),u);
    pq={FopConj,p,q};

end
