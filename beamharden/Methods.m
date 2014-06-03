classdef Methods < handle
    properties
        % for common purpose
        mainFunc = @(o) NCG_PR(o);
        alpha
        n = 3;
        fArray
        fVal
        coef=1;
        hArray
        difAlpha=1;
        difObj
        maxStepNum = 1e2;
        stepNum
        stepSize=0;
        converged = false;
        warned = false;
        % for NCG_PR

        % for SpaRSA
        Psi     % Psi and Psit are inv dwt transformations
        Psit
        sigma =0.01;    %
        t=-1;   % suppose be larger than L, largest eigenvalue of Hessian
        u       % coefficient of l1 norm term
        M = 5;  % keep the record of objective from the last M iterations
        oldCost
        cost;

        % for FISTA
        p = 0;
        ppp = 0;

    end
    methods
        function obj = Methods(n,alpha)
            if(nargin>0)
                obj.n = n;
            end
            obj.fArray = cell(obj.n,1);
            obj.fVal = zeros(obj.n,1);
            obj.hArray = cell(obj.n,1);
            obj.coef = zeros(obj.n,1);
            obj.alpha = alpha;
            obj.oldCost = zeros(obj.M,1);
        end

        function [f,g,h] = func(obj,aaa)
            if(nargout==3)
                g = 0;
                for i = 1:obj.n
                    [obj.fVal(i),gtemp,obj.hArray{i}] = obj.fArray{i}(aaa);
                    g = g+gtemp*obj.coef(i);
                end
                f = innerProd(obj.fVal(1:obj.n),obj.coef(1:obj.n));
                h = @(xxx,opt) hessian(xxx,opt);
            elseif(nargout==2)
                g = 0;
                for i = 1:obj.n
                    [obj.fVal(i),gtemp] = obj.fArray{i}(aaa);
                    g = g+gtemp*obj.coef(i);
                end
                f = innerProd(obj.fVal(1:obj.n),obj.coef(1:obj.n));
            else
                for i = 1:obj.n
                    obj.fVal(i) = obj.fArray{i}(aaa);
                end
                f = obj.fVal(1:obj.n)'*obj.coef(1:obj.n);
            end
            function hh = hessian(x,opt)
                hh=0;
                for j = 1:obj.n
                    hh = hh + obj.hArray{j}(x,opt)*obj.coef(j);
                end
                if(length(hh)==1 && hh<0)
                    warning(['\n%s: obj function is non-convex over alpha, '...
                        'x''*H*x=%g, replace it by 1'],class(obj),hh);
                    keyboard
                    obj.warned = true;
                    hh=1;
                end
            end
        end
        function set.M(obj,M)
            obj.M = M;
            obj.oldCost = zeros(obj.M,1);
        end
        function set.u(obj,u)
            obj.u = u;
            if(strcmpi(class(obj),'NCG_PR'))
                obj.coef(obj.n) = u;
            else
                obj.coef(obj.n+1)=u;
            end
        end
        function x= cg(c,hessianA,atHessianA,maxItr)
            % This function solve the problem 
            % min c'*x+1/2 atHessianA(x)
            % hessianA=hessian*x
            x=0; g=c; p=0; i=0;
            while(i<maxItr)
                i= i+1;
                preP=p; preG=g;
                g=c+hessianA(x);
                p=-g-g'*g/(preG'*preG)*preP;
                x=x-p*(p'*g/atHessianA(p));
            end
        end

        function oldJunk()
            if(t1==0 && p==1)
                if(interiorPointAlpha)
                    t1=min([1, abs(diff0'*difphi/pNorm(difphi)), abs(costA/costB)]);
                else t1=1;
                end
            end

            if(0)
                %afun=@(xxx,yyy) fhessianA(xxx);
                %[deltaAlpha,~]=bicg(afun,difAlpha,[],5);
                %afun=@(xxx,yyy) fhessianA(xxx);
                fhessianA=@(gAlpha) hessianA(gAlpha,weight,hphi*t1,Phi,Phit);
                fatHessianA=@(gAlpha) atHessianA(gAlpha,weight,t1*hphi,Phi,Phit);
                deltaAlpha = cg(-difAlpha,fhessianA,fatHessianA,1);
            end
            if(interiorPointAlpha)
                temp=find(deltaAlpha>0);
                if(isempty(temp)) maxStep=1;
                else maxStep=min(alpha(temp)./deltaAlpha(temp)); end
                maxStep=maxStep*0.99;
            end
            if(interiorPointAlpha)
                if(obj.deltaNormAlpha<1e-5)
                    if(t1 < 1e-10/length(alpha))
                        alphaReady=1;
                    else
                        t1=t1/10;
                        thresh1=thresh1/10;
                    end
                end
            end

        end
    end
end

