function proximal=sparseProximal(sparseType, prj_C, opt)
    % make an object to solve 0.5*||x-a||_2^2+u*TV(x)+I_C(x)
    % TV and C are specified via constructor see denoise for a and u

    switch(lower(sparseType))
        case 'iso'
            Psi=@(p) Psi_v(real(p))-Psi_h(imag(p));
            Psit=@(x) Psi_vt(x)-1i*Psi_ht(x);
            prj.op=@(p) p./max(1,abs(p));
            zeroInitSize=@(x) size(x);
            initStep='fixed';
            initStepSize=@(u) 1/16/u^2;
        case 'l1'
            Psi=@(p) Psi_v(p(1:end/2,:,:,:))+Psi_h(p(end/2+1:end,:,:,:));
            Psit=@(x) [Psi_vt(x); Psi_ht(x)];
            prj.op=@(p) min(max(p,-1),1);
            zeroInitSize=@(x) [size(x,1)*2, size(x,2)];
            initStep='fixed';
            initStepSize=@(u) 1/16/u^2;
        case 'realCustom'
            Psi=opt.Psi;
            Psit=opt.Psit;
            prj.op=@(p) min(max(p,-1),1);
            zeroInitSize=@(x) size(Psit(x));
            % Note that in special cases, such as DWT, initial step can
            % have better values to start from.
            initStep='bb';
            initStepSize=@(u) 1;
        otherwise
            error('error sparseType: %s\n',sparseType);
    end
    if(~exist('prj_C','var') || isempty(prj_C)) prj_C=@(x) x; end

    prj.penalty=@(x)0;
    prj.iterative=false;

    proximal.penalty = @penalty;
    proximal.op=@denoise;
    proximal.iterative=true;

    Psip=[]; x=[];

    function [x,itr,p]=denoise(a,u,thresh,maxItr,pInit)
        % Set default value for maxItr, thresh, and pInit, if needed.
        if(~exist('thresh','var')) thresh=1e-13; end
        if(~exist('maxItr','var')) maxItr=1e3; end
        if(~exist('pInit','var') || isempty(pInit)) pInit=zeroInitSize(a); end

        option.thresh=thresh; option.maxItr=maxItr;
        option.sol='PNPG'; option.adaptiveStep=true;
        option.debug=Debug(0);
        option.initStep='fixed';
        option.Lip=1/initStepSize(u);

        NLL=@(p) dualFunc(p,a,Psi,Psit,u,prj_C);

        out=pnpg(NLL,prj,pInit,option);

        p=out.x;
        Psip=Psi(p);
        itr=out.itr;
        x=prj_C(a-u*Psip);
    end
    function pen=penalty(z)
        if(nargin==0)
            % Psi(p) is real
            pen=innerProd(x,Psip);
        else
            pen=pNorm(Psit(z),1);
        end
    end
end

function [f,g,h] = dualFunc(p,a,Psi,Psit,u,prj_C)
    % minimize 0.5*||a-u*real(Psi(p))||_2^2-0.5*||X-a+u*real(Psi(p))||_2^2
    % subject to ||p||_infty <= 1
    % where X=prj_C( a-u*real(Psi(p)) ), p and Psi may be complex.
    Qp=a-u*Psi(p); % is real
    x=prj_C(Qp);   % is real
    f=(sqrNorm(Qp)-sqrNorm(x-Qp))/2;
    if(nargout>1)
        g=-u*Psit(x);
        if(nargout>2)
            h=[];
        end
    end
end

% if edit the following, update TV.[A,B]t\=
function x = Psi_v(p)
    [I,~]=size(p);
    p(I,:)=0;
    x=[p(1,:); p(2:I,:)-p(1:I-1,:)];
end
function p = Psi_vt(x)
    [I,J]=size(x);
    p=[x(1:I-1,:)-x(2:I,:);zeros(1,J)];
end
function x = Psi_h(q)
    [~,J]=size(q);
    q(:,J)=0;
    x=[q(:,1), q(:,2:J)-q(:,1:J-1)];
end
function q = Psi_ht(x)
    [I,J]=size(x);
    q=[x(:,1:J-1)-x(:,2:J), zeros(I,1)];
end

