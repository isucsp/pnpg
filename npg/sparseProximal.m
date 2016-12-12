function proximal=sparseProximal(sparseType, prj_C, opt, option)
    % make an object to solve 0.5*||x-a||_2^2+u*TV(x)+I_C(x)
    % TV and C are specified via constructor see denoise for a and u

    switch(lower(sparseType))
        case 'iso'
            Psi=@(p) Psi_v(real(p))-Psi_h(imag(p));
            Psit=@(x) Psi_vt(x)-1i*Psi_ht(x);
            prj.op=@(p) p./sqrt(max(1,real(p).^2+imag(p).^2));
            zeroInit=@(x) zeros(size(x));
            initStep='fixed';
            Lipschitz=@(u) 8*u^2;
        case 'l1'
            Psi=@(p) Psi_v(p(1:end/2,:,:,:))+Psi_h(p(end/2+1:end,:,:,:));
            Psit=@(x) [Psi_vt(x); Psi_ht(x)];
            prj.op=@(p) min(max(p,-1),1);
            zeroInit=@(x) zeros([size(x,1)*2, size(x,2)]);
            initStep='fixed';
            Lipschitz=@(u) 8*u^2;
        case 'realCustom'
            Psi=opt.Psi;
            Psit=opt.Psit;
            prj.op=@(p) min(max(p,-1),1);
            zeroInit=@(x) zeros(size(Psit(x)));
            % Note that in special cases, such as DWT, initial step can
            % have better values to start from.
            initStep='bb';
            Lipschitz=@(u) 1;
        otherwise
            error('error sparseType: %s\n',sparseType);
    end
    if(~exist('prj_C','var') || isempty(prj_C)) prj_C=@(x) x; end

    prj.penalty=@(x)0;
    prj.iterative=false;

    proximal.penalty = @penalty;
    proximal.op=@denoise;
    proximal.iterative=true;

    if(~exist('option','var') || ~isfield(option,'adaptiveStep')) option.adaptiveStep=true; end
    if(~isfield(option,'debugLevel')) option.debugLevel=0; end
    if(~isfield(option,'debugLevel')) option.outLevel=0; end
    if(~isfield(option,'initStep')) option.initStep=initStep; end
    if(~isfield(option,'initStep')) option.minItr=1; end

    function [x,itr,p,out]=denoise(a,u,thresh,maxItr,pInit)
        % Set default value for maxItr, thresh, and pInit, if needed.
        if(~exist('thresh','var')) thresh=1e-13; end
        if(~exist('maxItr','var')) maxItr=1e3; end
        if(~exist('pInit','var') || isempty(pInit)) pInit=zeroInit(a); end

        option.thresh=thresh; option.maxItr=maxItr;
        option.Lip=Lipschitz(u);

        NLL=@(p) dualFunc(p,a,Psi,Psit,u,prj_C);

        out=pnpg(NLL,prj,pInit,option);

        itr=out.itr;
        p=out.x;
        x=prj_C(a-u*Psi(p));
    end
    function pen=penalty(z)
        pen=pNorm(Psit(z),1);
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

% If edit the following, update TV.[A,B]t\=
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

