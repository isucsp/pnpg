function [sp,penalty]=sparseProximal(sparseType, prj_C, opt)
    % make an object to solve 0.5*||x-a||_2^2+u*TV(x)+I_C(x)
    % TV and C are specified via constructor see denoise for a and u

    switch(lower(sparseType))
        case 'iso'
            Psi=@(p) Psi_v(real(p))-Psi_h(imag(p));
            Psit=@(x) Psi_vt(x)-1i*Psi_ht(x);
            opt.proximalOp=@(p) p./max(1,abs(p));
            zeroInitSize=@(x) size(x);
            initStepSize=@(u) 1/16/u^2;
        case 'l1'
            Psi=@(p) Psi_v(p(1:end/2,:,:,:))+Psi_h(p(end/2+1:end,:,:,:));
            Psit=@(x) [Psi_vt(x); Psi_ht(x)];
            opt.proximalOp=@(p) min(max(p,-1),1);
            zeroInitSize=@(x) [size(x,1)*2, size(x,2)];
            initStepSize=@(u) 1/16/u^2;
        case 'realCustom'
            Psi=opt.Psi;
            Psit=opt.Psit;
            opt.proximalOp=@(p) min(max(p,-1),1);
            zeroInitSize=@(x) size(Psit(x));
            initStepSize=@(u) 1;
        otherwise
            error('error sparseType: %s\n',sparseType);
        end
        penalty = @(x) pNorm(Psit(x),1);
    end

    if(exist('prj_C','var')) prj_C=prj_C; else prj_C=@(x) x; end

    sp=@denoise;

    function [x,itr,p]=denoise(a,u,thresh,maxItr,p)
        % Set default value for maxItr, thresh, and p, if needed.
        if(~exist('thresh','var')) thresh=1e-13; end
        if(~exist('maxItr','var')) maxItr=1e3; end
        if(~exist('p','var') || isempty(p)) p=zeroInitSize(a); end

        opt.thresh=thresh; opt.maxItr=maxItr;

        opt.NLL=@(p) dualTvFunc(p,a,Psi,Psit,u,prj_C);

        ??? opt.continuation=false; opt.sol='PNPG'; opt.adaptiveStep=true;
        debug=Debug(0);

        initStepSize(u);

        %call npg here
        % No need to set u
        %sol=NPG(1,init,maxItr,[],proximalOp);
        % ??? out=solver([],[],[],[],[],init,opt);

        p=out.alpha;
        x=prj_C(a-u*Psi(p));
        itr=out.p;
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
end

% if edit the following, update TV.[A,B]t\=
function x = Psi_v(p)
    [I,J]=size(p);
    p(I,:)=0;
    x=[p(1,:); p(2:I,:)-p(1:I-1,:)];
end
function p = Psi_vt(x)
    [I,J]=size(x);
    p=[x(1:I-1,:)-x(2:I,:);zeros(1,J)];
end
function x = Psi_h(q)
    [I,J]=size(q);
    q(:,J)=0;
    x=[q(:,1), q(:,2:J)-q(:,1:J-1)];
end
function q = Psi_ht(x)
    [I,J]=size(x);
    q=[x(:,1:J-1)-x(:,2:J), zeros(I,1)];
end

