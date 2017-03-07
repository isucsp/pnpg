function u=uBound(Psi,Psit,tvType,prj_NC,xstar,g,opt)
%
% Upper-Bounding the Regularization Constant for the following problem
%
%               minimize     L(x) + u*||\Psi^H x||_1                   (1)
%               subject to   x \in C
%
% Phi and Phit: function handle for Phi and Phi^H
% tvType:       can be [], 'wav', 'iso', and 'ani'
% prj_NC:       function handle for the projection to the normal cone of
%               convex set C at x_0, where x_0 is the solution of (1) when
%               u -> \infty.
% xstar:        equivalent to x_0.
% g:            the gradient of L(x) at xstar, i.e. x_0
% opt:          options
%
% Example: 
%
%       For the problem L(x)=0.5*||x||_2^2 with C={x| ||x-[2;0]||_2^2<=2 },
%       suppose Psi and Psit are operators for 1-d TV of a signal x, such
%       that size(x)=[2,1].  The following code returns the upper bound of
%       u:
%      
%               xstar=ones(2,1);
%               Psi=@(x)[x;-x];
%               Psit=@(x) x(1)-x(2);
%               %Ncx=[-a,a] with a>=0
%               Pncx=@(x) min(0,(x(1)-x(2))/2)*[1;-1];
%               g=[1;1];
%               uTrue=Inf
%               u=uBound(Psi,Psit,[],Pncx,xstar,g);
%
%
% See Also:     uBoundEx/{exa1,exa2,exa3,slGaussEx,tv_Bound}
% Author:       Renliang Gu
% Reference: 
%       R. Gu and A. Dogandžić. Upper-Bounding the Regularization Constant
%       for Convex Sparse Signal Reconstruction.


% Theoretically, the condition is 
%         0 \in g + N_C(x*)
% where N_C is the normal cone of C at xstar.
if(norm(reshape(g+prj_NC(-g),[],1),2)<eps)
    u=0;
    return;
end

if(~exist('opt','var') || ~isfield(opt,'maxItr')) opt.maxItr=5e3; end

% it is sure that normG~=0
normG=norm(g(:),2);
g=g/normG;

v=1;
vg=v*g;
t=zeros(size(xstar));
z=zeros(size(t));
if(isempty(tvType))
    w=zeros(size(Psit(xstar)));
    realPsi=Psi;
    realPsit=Psit;
else
    switch(lower(tvType))
        case {'iso','l1'}
            w=zeros(size([TV.Phi_vt(xstar); TV.Phi_ht(xstar)]));
            realPsi=@(w) TV.Phi_v(w(1:end/2,:))+TV.Phi_h(w(end/2+1:end,:));
            realPsit=@(x) [TV.Phi_vt(x); TV.Phi_ht(x)];
        case {'wav'}
            w=zeros(size(Psit(xstar)));
            realPsi=Psi;
            realPsit=Psit;
    end
end

rho=1;

incCnt=0;
decCnt=0;
ii=0;
EndCnt=0;

vgtz=zeros(size(g));

strlen=0;

PsitXstar=realPsit(xstar);
if(isempty(tvType))
    lb=-ones(size(w)); ub=ones(size(w));
    ww=sign(PsitXstar); A=(ww~=0);
    lb(A)=ww(A); ub(A)=ww(A);
    proj=@(w) max(min(w,ub),lb);
else
    switch(lower(tvType))
        case {'iso'}
            proj=Prj_G(PsitXstar);
        case {'l1','wav'}
            lb=-ones(size(w)); ub=ones(size(w));
            ww=sign(PsitXstar); A=(ww~=0);
            lb(A)=ww(A); ub(A)=ww(A);
            proj=@(w) max(min(w,ub),lb);
    end
end

% Lip = 8 for 2-d TV
% Lip = 4 for 1-d TV
% Lip = 1 for wavelet
if(isempty(tvType))
    opts.printEvery=inf;
    opts.maxTotalIts=5e3;
    opts.maxIts=1000;
    opts.factr=0; %1e-6/eps;
else
    switch(lower(tvType))
        case {'iso','l1'}
            if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=false; end
            if(~isfield(opt,'backtracking')) opt.backtracking=false; end
            if(~isfield(opt,'debugLevel')) opt.debugLevel=0; end
            if(~isfield(opt,'outLevel')) opt.outLevel=0; end
            if(~isfield(opt,'Lip'))
                opt.Lip=8;
                if(min(size(g))==1) opt.Lip=4; end
            end
            if(~isfield(opt,'initStep')) opt.initStep='fixed'; end
            proximal.exact=true; proximal.val=@(x)0; proximal.prox=proj;

            opts.printEvery=inf;
            opts.maxTotalIts=5e3;
            opts.maxIts=1000;
            opts.factr=0; %1e-6/eps;
        case {'wav'}
            proximal.exact=true; proximal.val=@(x)0; proximal.prox=proj;
    end
end

fprintf('\n%s\n', repmat( '=', 1, 80 ) );
str=sprintf('ADMM for uBound, Type: ');
if(isempty(tvType))
    str=[str 'Sparsifying Matrix Psi'];
else
    switch(lower(tvType))
        case 'iso'
            str=[str 'ISO_TV'];
        case 'l1'
            str=[str 'L1_TV'];
        case 'wav'
            str=[str 'Wavelet Psi'];
    end
end

fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
fprintf('%s\n', repmat('=',1,80));
str=sprintf( '%5s','Itr');
str=sprintf([str ' %12s'],'u');
str=sprintf([str ' %12s'], 'PrimRes');
str=sprintf([str ' %12s'], 'DualRes1');
str=sprintf([str ' %10s'], 'DualRes2');
str=sprintf([str ' %12s'], 'DualGap');
str=sprintf([str ' %4s'], 'iItr');
str=sprintf([str ' %4s'], 'rho');
fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );

cntThresh=100;
while(EndCnt<3 && ii<=opt.maxItr)

    ii=ii+1;

    preV=v; preT=t; preVgtz=vgtz;

    vgtz=(vg+t+z);
    NLL=@(x) Utils.linearModel(x,realPsi,realPsit,-vgtz);
    if isempty(tvType)
        % objective: 0.5*||Psi(w)+t+z+v*g||_F^2
        % subject to the [-1,1] box
        opts.pgtol=max(1e-3*relativeDif(preVgtz,vgtz),1e-14);
        opts.x0=w;
        [w,~,info]=lbfgsb(NLL,lb,ub,opts);
        innerItr=info.iterations;
        %numItr=info.iterations; %disp(info); strlen=0;
    elseif strcmpi(tvType,'iso') || strcmpi(tvType,'l1')
        %opt.debugLevel=2; opt.outLevel=1;
        opt.thresh=max(1e-4*relativeDif(preVgtz,vgtz),1e-13);
        out = pnpg(NLL,proximal,w,opt);
        w=out.x; innerItr=out.itr;
    elseif strcmpi(tvType,'wav')
        w=proximal.prox(-realPsit(vgtz));
        innerItr=-1;
    end

    realPsi_w=realPsi(w);

    % Since ||g||=1, v=(1-rho*g'*(realPsi_w+t+z)) / (rho*(g'*g))
    % reduces to
    v=1/rho-sum(reshape(g.*(realPsi_w+t+z),[],1)); vg=v*g;
    t=prj_NC( -vg-realPsi_w-z );

    z=z+realPsi_w+vg+t;

    pRes=norm(reshape(vg+realPsi_w+t,[],1),2);
    dRes1=abs(preV-v);  %norm(g*(preV-v));
    dRes2=norm(preT-t);
    dRes=max(dRes1,dRes2);
    gap=sum(reshape(z.*(vg+realPsi_w+t),[],1));

    str=sprintf('%5d %12g %12g %12g %10g %12g %4g %4g',ii, normG/v, pRes,...
        dRes1, dRes2, gap,innerItr, rho);
    if(~isfield(opt,'debugLevel') || opt.debugLevel==0)
        if(strlen==0 || (mod(ii-1,100)==0 || (ii<=100 && mod(ii-1,10)==0) || ii-1<10))
            fprintf('\n%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
    end
    strlen = length(str);

    if(dRes>10*pRes) decCnt = decCnt+1; else decCnt = 0; end
    if(dRes<pRes/10) incCnt = incCnt+1; else incCnt = 0; end

    if(decCnt>cntThresh)
        rho=rho/2; z=z*2; decCnt=0; cntThresh=cntThresh*1.6;
    end
    if(incCnt>cntThresh)
        rho=rho*2; z=z/2; incCnt=0; cntThresh=cntThresh*1.6;
    end

    if(max(pRes,dRes)<1e-9)
        EndCnt=EndCnt+1;
    else
        EndCnt=0;
    end
end
u=normG/v;
fprintf('\n\n');

function y = Prj_G(PsitXstar)
    [I,~]=size(PsitXstar); I=I/2;

    %mag=sqrt(max(1,[p.^2;zeros(1,J)]+[q.^2, zeros(I,1)]));
    pp=PsitXstar(1:I,:); qq=PsitXstar(I+1:end,:);
    mag=sqrt(pp.^2+qq.^2); nonzeroMag=(mag~=0);
    pp(nonzeroMag)=pp(nonzeroMag)./mag(nonzeroMag);
    qq(nonzeroMag)=qq(nonzeroMag)./mag(nonzeroMag);
    
    y=@(w) prj(w);

    function y = prj(w)

        p=w(1:I,:); q=w(I+1:end,:);
        mag2=max(1,sqrt(p.^2+q.^2));
        p=p./mag2; q=q./mag2;
        p(nonzeroMag)=pp(nonzeroMag);
        q(nonzeroMag)=qq(nonzeroMag);

        y=[p; q];
    end
end
end

