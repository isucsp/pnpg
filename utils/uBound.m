function u=uBound(Psi,Psit,prj_NC,xstar,g)
% 
% g is the gradient of L(x) at xstar
% author: Renliang Gu
%

% Theoretically, the condition is 
%         0 \in g + N_C(x*)
% where N_C is the normal cone of C at xstar.
if(norm(reshape(g+prj_NC(-g),[],1),2)<eps)
    u=0;
    return;
end

% it is sure that normG~=0
normG=norm(g(:),2);
g=g/normG;

v=1;
vg=v*g;
t=zeros(size(xstar));
if(isstruct(Psit))
    w=zeros(size([Psit.r(xstar); -Psit.i(xstar)]));
    realPsi=@(w) Psi.r(w(1:end/2,:))-Psi.i(w(end/2+1:end,:));
    realPsit=@(x) [Psit.r(x); -Psit.i(x)];
else
    w=zeros(size(Psit(xstar)));
    realPsi=Psi;
    realPsit=Psit;
end
z=zeros(size(t));
realPsi_w=realPsi(w);

rho=1;

cnt=0;
ii=0;
EndCnt=0;

dRes=1e-4;

strlen=0;

PsitXstar=realPsit(xstar);
if(isstruct(Psi))
    opts.printEvery=inf;
    opts.maxTotalIts=5e3;
    opts.maxIts=1000;
    opts.factr=0; %1e-6/eps;

    lb=-ones(size(w)); ub=ones(size(w));
    ww=sign(PsitXstar); A=(ww~=0);
    lb(A)=ww(A); ub(A)=ww(A);
    proj=@(w) max(min(w,ub),lb);
else
    tt=PsitXstar(1:end/2,:)+1i*PsitXstar(end/2+1:end,:);
    zeroIdx=(tt==0);
    PsitXstarNonZeroNormalized=conj(tt(~zeroIdx))./abs(tt(~zeroIdx));
    proj=@(w) projection(w,zeroIdx,PsitXstarNonZeroNormalized);
end

while(EndCnt<3 && ii<=5e3)

    ii=ii+1; cnt=cnt+1;

    preV=v; preT=t; preW=w; preZ=z; preRealPsi_w=realPsi_w;

    if(isstruct(Psi))
        % objective: 0.5*||Psi(w)+t+z+v*g||_F^2
        % subject to the [-1,1] box
        opts.pgtol=max(1e-2*dRes,1e-14);
        opts.x0=w;
        [w,cost,info]=lbfgsb(@(x) quadBox(x,vg+t+z,realPsi,realPsit),...
            lb,ub,opts);
        %numItr=info.iterations; %disp(info); strlen=0;
    else
        if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=false; end
        if(~isfield(opt,'backtracking')) opt.backtracking=false; end
        if(~isfield(opt,'debugLevel')) opt.debugLevel=0; end
        if(~isfield(opt,'outLevel')) opt.outLevel=0; end

        % Lip = 8 for 2-d TV
        % Lip = 4 for 1-d TV
        % Lip = 1 for wavelet

        vgtz=(vg+t+z); NLL=@(x) Utils.linearModel(x,realPsi,realPsit,-vgtz);
        proximal.exact=true; proximal.val=@(x)0; proximal.prox=proj;
        out = pnpg(NLL,proximal,w,opt)
        w=out.x;
    end
    realPsi_w=realPsi(w);

    % Since ||g||=1, v=(1-rho*g'*(realPsi_w+t+z)) / (rho*(g'*g))
    % reduces to
    v=1/rho-g'*(realPsi_w+t+z); vg=v*g;
    t=prj_NC( -vg-realPsi_w-z );

    z=z+realPsi_w+vg+t;

    pRes=norm(vg+realPsi_w+t,2);
    %dRes1=g'*(preRealPsi_w-realPsi_w+preT-t);
    %dRes2=norm(preRealPsi_w-realPsi_w,2)^2;
    dRes1=abs(preV-v);  %norm(g*(preV-v));
    dRes2=norm(preT-t);
    dRes=[max(dRes1,dRes2)];
    gap=z'*(vg+realPsi_w+t);

    str=sprintf('itr=%d, u=%g pRes=%g dRes1=%g dRes2=%g gap=%g rho=%g                 ',ii, normG/v, pRes,...
        dRes1, dRes2, gap, rho);
    if(strlen==0 || (mod(ii-1,100)==0 || (ii<=100 && mod(ii-1,10)==0) || ii-1<10))
        fprintf('\n%s',str);
    else
        fprintf([repmat('\b',1,strlen) '%s'],str);
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
function y = projection(x,zeroIdx,PsitXstarNonZeroNormalized)
    x=x(1:end/2,:)+1i*x(end/2+1:end,:);
    x(zeroIdx)=x(zeroIdx)./max(1,abs(x(zeroIdx)));
    x(~zeroIdx)=(PsitXstarNonZeroNormalized);
    y=[real(x); imag(x)];
end
end

