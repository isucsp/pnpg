function u=uBound(Psi,Psit,tvType,prj_NC,xstar,g,opt)
% 
% g is the gradient of L(x) at xstar
%
%
%
% Examples: 
%
% author: Renliang Gu
%

% Theoretically, the condition is 
%         0 \in g + N_C(x*)
% where N_C is the normal cone of C at xstar.
if(norm(reshape(g+prj_NC(-g),[],1),2)<eps)
    u=0;
    return;
end

if(~exist('opt','var')) opt=[]; end

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
    w=zeros(size([TV.Phi_vt(xstar); TV.Phi_ht(xstar)]));
    realPsi=@(w) TV.Phi_v(w(1:end/2,:))+TV.Phi_h(w(end/2+1:end,:));
    realPsit=@(x) [TV.Phi_vt(x); TV.Phi_ht(x)];
end
realPsi_w=realPsi(w);

rho=1;

cnt=0;
ii=0;
EndCnt=0;

dRes=1e-4;

strlen=0;

PsitXstar=realPsit(xstar);
switch(lower(tvType))
    case {'iso','l1'}
        proj=Prj_G(PsitXstar);
    otherwise
        lb=-ones(size(w)); ub=ones(size(w));
        ww=sign(PsitXstar); A=(ww~=0);
        lb(A)=ww(A); ub(A)=ww(A);
        proj=@(w) max(min(w,ub),lb);
end

% Lip = 8 for 2-d TV
% Lip = 4 for 1-d TV
% Lip = 1 for wavelet
switch(lower(tvType))
    case {'iso','l1'}
        if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=false; end
        if(~isfield(opt,'backtracking')) opt.backtracking=false; end
        if(~isfield(opt,'debugLevel')) opt.debugLevel=0; end
        if(~isfield(opt,'outLevel')) opt.outLevel=0; end
        if(~isfield(opt,'Lip')) opt.Lip=8; end
        if(~isfield(opt,'initStep')) opt.initStep='fixed'; end
        proximal.exact=true; proximal.val=@(x)0; proximal.prox=proj;
    otherwise
        opts.printEvery=inf;
        opts.maxTotalIts=5e3;
        opts.maxIts=1000;
        opts.factr=0; %1e-6/eps;
end

fprintf('\n%s\n', repmat( '=', 1, 80 ) );
str=sprintf('ADMM for uBound, Type: ');
switch(lower(tvType))
    case 'iso'
        str=[str 'ISO_TV'];
    case 'l1'
        str=[str 'L1_TV'];
    otherwise
        str=[str 'Sparsifying Matrix Psi'];
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

while(EndCnt<3 && ii<=5e3)

    ii=ii+1; cnt=cnt+1;

    preV=v; preT=t; preW=w; preZ=z; preRealPsi_w=realPsi_w;

    switch(lower(tvType))
        case {'iso','l1'}
            vgtz=(vg+t+z); NLL=@(x) Utils.linearModel(x,realPsi,realPsit,-vgtz);
            %opt.debugLevel=2; opt.outLevel=1;
            out = pnpg(NLL,proximal,w,opt);
            w=out.x;
            innerItr=out.itr;
            %opt.thresh=1e-3*relativeDif(preW,w);
        otherwise
            % objective: 0.5*||Psi(w)+t+z+v*g||_F^2
            % subject to the [-1,1] box
            opts.pgtol=max(1e-2*dRes,1e-14);
            opts.x0=w;
            [w,cost,info]=lbfgsb(@(x) quadBox(x,vg+t+z,realPsi,realPsit),...
                lb,ub,opts);
            innerItr=info.iterations;
            %numItr=info.iterations; %disp(info); strlen=0;
    end

    realPsi_w=realPsi(w);

    % Since ||g||=1, v=(1-rho*g'*(realPsi_w+t+z)) / (rho*(g'*g))
    % reduces to
    v=1/rho-g'*(realPsi_w+t+z); vg=v*g;
    t=prj_NC( -vg-realPsi_w-z );

    z=z+realPsi_w+vg+t;

    pRes=norm(reshape(vg+realPsi_w+t,[],1),2);
    %dRes1=g'*(preRealPsi_w-realPsi_w+preT-t);
    %dRes2=norm(preRealPsi_w-realPsi_w,2)^2;
    dRes1=abs(preV-v);  %norm(g*(preV-v));
    dRes2=norm(preT-t);
    dRes=[max(dRes1,dRes2)];
    gap=z'*(vg+realPsi_w+t);

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

    if(cnt>100) % prevent excessive back and forth adjusting
        if(dRes>10*pRes)
            rho=rho/2; z=z*2; cnt=0;
        elseif(dRes<pRes/10)
            rho=rho*2; z=z/2; cnt=0;
        end
    end

    if(max(pRes,dRes)<1e-9)
        EndCnt=EndCnt+1;
    else
        EndCnt=0;
    end
end
u=normG/v;
fprintf('\n\n');

function [f,g] = quadBox(x,y,Psi,Psit)
    r=Psi(x)+y;
    f=0.5*norm(r,2)^2;
    if(nargout>1)
        g=Psit(r);
    end
end
function [f,g] = compBox(x,y,Psi,Psit)
    r=Psi.r(x(:,1))+Psi.i(x(:,2))+y;
    f=0.5*sqrNorm(r);
    if(nargout>1)
        g=[Psit.r(r), Psit.i(r)];
    end
end
function y = Prj_G(PsitXstar)
    [I,J]=size(PsitXstar); I=I/2;

    %mag=sqrt(max(1,[p.^2;zeros(1,J)]+[q.^2, zeros(I,1)]));
    pp=PsitXstar(1:I,:); qq=PsitXstar(I+1:end,:);
    mag=sqrt(pp.^2+qq.^2); nonzeroMag=(mag~=0);
    pp(nonzeroMag)=pp(nonzeroMag)/mag;
    qq(nonzeroMag)=qq(nonzeroMag)/mag;
    
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

