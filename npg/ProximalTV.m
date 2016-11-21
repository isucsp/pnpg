classdef ProximalTV < Proximal
    properties
        opt
        prj_C=@(x)x;
        Psi
        Psit
        out % put this property to be private later
    end
    properties (Access=private)
        zeroInitSize
        temp
    end
    properties (Dependent)
        tvType
    end
    methods
        % make an object to solve 0.5*||x-a||_2^2+u*TV(x)+I_C(x)
        % TV and C are specified via constructor see denoise for a and u
        function obj = ProximalTV(tvType,prj_C)
            obj.tvType=tvType;
            if(exist('prj_C','var'))
                obj.prj_C=prj_C;
            else
                obj.prj_C=@(x) x;
            end
            obj.iterative=true;
            obj.preConfig();
        end
        % don't forget to set default value for maxItr, thresh
        function set.tvType(obj,tvType)
            switch(lower(tvType))
                case 'iso'
                    obj.Psi=@(p) TV.A(real(p))-TV.B(imag(p));
                    obj.Psit=@(x) TV.At(x)-1i*TV.Bt(x);
                    obj.opt.proximalOp=@(p) p./max(1,abs(p));
                    obj.zeroInitSize=@(x) size(x);
                case 'l1'
                    obj.Psi=@(p) TV.A(p(1:end/2,:,:,:))+TV.B(p(end/2+1:end,:,:,:));
                    obj.Psit=@(x) [TV.At(x); TV.Bt(x)];
                    obj.opt.proximalOp=@(p) min(max(p,-1),1);
                    obj.zeroInitSize=@(x) [size(x,1)*2, size(x,2)];
                otherwise
                    error('error tvType: %s\n',tvType);
            end
            obj.regFunc = @(X) tlv(X,tvType);
        end
        function setInit(obj)
            obj.init=obj.out.alpha;
        end
        function x=denoise(obj,a,u)
            obj.config(a,u);
            obj.out=solver([],[],[],[],[],obj.init,obj.opt);
            obj.temp=obj.Psi(obj.out.alpha);
            obj.x=obj.prj_C(a-u*obj.temp);
            obj.steps=obj.out.p;
            x=obj.x;
        end
        function y=getPenalty(obj,X)
            if(~exist('X','var'))
                y=innerProd(obj.x,obj.temp);
            else
                y=pNorm(obj.Psit(X),1);
            end
        end
        %     ??? check complex variable of solver
        function config(obj,a,u)
            obj.opt.thresh=obj.thresh; obj.opt.maxItr=obj.maxItr;
            if(any(size(obj.init)-obj.zeroInitSize(a)~=0))
                obj.init=zeros(obj.zeroInitSize(a));
            end
            obj.opt.NLL=@(p) dualTvFunc(p,a,obj.Psi,obj.Psit,u,obj.prj_C);
            function [f,g,h] = dualTvFunc(p,a,Psi,Psit,u,prj_C)
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
        function preConfig(obj)
            obj.opt.continuation=false; obj.opt.alphaStep='PNPG'; obj.opt.adaptiveStep=true;
        end
    end
end
