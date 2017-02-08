classdef TV < handle
methods(Static,Access=private)
function out = weight(i,j,k,l)
    %               /k+l-2\ /i+j-k-l-1\
    %               \ k-1 / \   j-l   /
    % approximate   -------------------
    %                     /i+j-2\
    %                     \ i-1 /
    % using Stirling's approximation

    if(i<k || j<l) 
        global strlen
        strlen=0;
        fprintf('TV.weight(i=%d, j=%d, k=%d, l=%d) needs larger i and j\n',i,j,k,l);
        return;
    end
    out = -0.5*log(2*pi);
    t=k+l-2;     out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=k-1;       out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=l-1;       out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);

    t=i+j-k-l-1; out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=i-1;       out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=j-1;       out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=i-k-1;     out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=j-l;       out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    t=i+j-2;     out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
    out = exp(out); 
end
function out = U1(g)
    [I,J]=size(g);
    G = zeros(2*(I+1),2*(J+1));
    G(2:I+1,2:J+1)=g;

    % size of G is 2(I+1) x 2(J+1)
    % index of G is from 0 to 2(I+1)-1,   0 to 2(J+1)-1

    p1p=zeros(2*(I+1),2*(J+1));
    p1p(1:I,1)=1;
    P1=ifft2(fft2(G).*fft2(p1p));
    P1=P1(2:I+1,2:J+1);

    s=sum(P1,2);
    H = zeros(I,J);
    H(:,1:J-1)=H(:,1:J-1)-repmat(s,1,J-1);
    H(:,2:J  )=H(:,2:J  )-repmat(s,1,J-1);
    H=H-2*(I-1)*P1;
    for s=1:I
        H(s,:)=(2*s-1)*P1(I,:);
    end
    H=H/(I+J-2)/2;
    out=max(abs(H(:)));
end

end % methods(Static,Access=private)

methods(Static)
function [newX,innerSearch]=denoise(x,u,innerThresh,maxInnerItr,maskmt,tvType,lb,ub,nRow,nCol)
    if(~exist('tvType','var')) tvType='iso'; end
    if(~exist('lb','var')) lb=0; end
    if(~exist('ub','var')) ub=inf; end
    if(~exist('nRow','var')) nRow=sqrt(length(x(:))); end
    if(~exist('nCol','var')) nCol=length(x(:))/nRow; end
    pars.print = 0;
    pars.tv = tvType;
    pars.MAXITER = maxInnerItr;
    pars.epsilon = innerThresh; 

    if(~exist('maskmt','var') || isempty(maskmt))
        mask.a=@(xx) xx(:);
        mask.b=@(xx) reshape(xx,nRow,nCol);
    else
        maskIdx=find(maskmt~=0);
        n=size(maskmt);
        mask.a=@(xx) maskFunc(xx,maskIdx);
        mask.b=@(xx) maskFunc(xx,maskIdx,n);
    end

    [newX,innerSearch]=denoise_bound_mod(mask.b(x),u,lb,ub,pars);
    newX=mask.a(newX);
end

% if edit the following, update sparseProximal
function x = A(p)
    [I,J]=size(p);
    p(I,:)=0;
    x=[p(1,:); p(2:I,:)-p(1:I-1,:)];
end
function p = At(x)
    [I,J]=size(x);
    p=[x(1:I-1,:)-x(2:I,:);zeros(1,J)];
end
function x = B(q)
    [I,J]=size(q);
    q(:,J)=0;
    x=[q(:,1), q(:,2:J)-q(:,1:J-1)];
end
function q = Bt(x)
    [I,J]=size(x);
    q=[x(:,1:J-1)-x(:,2:J), zeros(I,1)];
end
function invP = inv_AtA(p)
    [I,J]=size(p); I=I+1;
    [k,ii]=meshgrid(1:I-1);
    matrix=min(k,ii).*min(I-k,I-ii);
    invP=matrix*p/I;
end
function invQ = inv_BtB(q)
    [I,J]=size(q); J=J+1;
    [k,j]=meshgrid(1:J-1);
    matrix=min(k,j).*min(J-k,J-j);
    invQ=q*matrix/J;
end

function [p,q]=isoPrj(p,q)
    [I,J]=size(p);
    I=I+1;
    %mag=sqrt(max(1,[p.^2;zeros(1,J)]+[q.^2, zeros(I,1)]));
    mag=sqrt(max(1,p(:,1:J-1).^2+q(1:I-1,:).^2));
    p(:,1:J-1)=p(:,1:J-1)./mag; p(:,J)=min(max(p(:,J),-1),1);
    q(1:I-1,:)=q(1:I-1,:)./mag; q(I,:)=min(max(q(I,:),-1),1);
end
function [p,q]=l1Prj(p,q)
    p=min(max(p,-1),1);
    q=min(max(q,-1),1);
end

function x = Psi(p,q)
    [I,J]=size(p); I=I+1;
    x=diff([zeros(1,J); p; zeros(1,J)],1,1) ...
        + diff([zeros(I,1), q, zeros(I,1)],1,2);
end

function [p,q] = Psit(x)
    p=-diff(x,1,1);
    q=-diff(x,1,2);
end

function x = Psi_v(p)
    %[I,J]=size(p);
    %I=I+1;
    %x=[p; zeros(1,J)];
    %x(2:I,:)=x(2:I,:)-p(1:I-1,:);
    %x=[p(1,:); diff(p,1,1); -p(end,:)];
    J=size(p,2);
    x=diff([zeros(1,J); p; zeros(1,J)],1,1);
end
function x = Phi_v(p)
    [I,J]=size(p);
    p(I,:)=0;
    %x=[p; zeros(1,J)];
    %x(2:I,:)=x(2:I,:)-p(1:I-1,:);
    %x=[p(1,:); diff(p,1,1); -p(end,:)];
    x=diff([zeros(1,J); p],1,1);
end
function p = Psi_vt(x)
    %[I,~]=size(x);
    %p=x(1:I-1,:)-x(2:I,:);
    p=-diff(x,1,1);
end
function p = Phi_vt(x)
    [I,J]=size(x);
    %p=x(1:I-1,:)-x(2:I,:);
    p=[-diff(x,1,1); zeros(1,J)];
end

function x = Psi_h(q)
    %[I,J]=size(q);
    %J=J+1;
    %x=[q zeros(I,1)];
    %x(:,2:J)=x(:,2:J)-x(:,1:J-1);
    I=size(q,1);
    x=diff([zeros(I,1), q, zeros(I,1)],1,2);
end
function x = Phi_h(q)
    [I,J]=size(q);
    q(:,J)=0;
    %x=[q zeros(I,1)];
    %x(:,2:J)=x(:,2:J)-x(:,1:J-1);
    x=diff([zeros(I,1), q],1,2);
end
function q = Psi_ht(x)
    %[~,J]=size(x);
    %q=x(:,1:J-1)-x(:,2:J);
    q=-diff(x,1,2);
end
function q = Phi_ht(x)
    [I,~]=size(x);
    %q=x(:,1:J-1)-x(:,2:J);
    q=[-diff(x,1,2), zeros(I,1)];
end

end % methods(Static)
end % classdef TV < handle

