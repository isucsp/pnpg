classdef ConjGrad < handle
    properties
        n = 3;
        cgType = 'pr';
        alpha
        fArray
        fVal
        hArray
        coef
        thresh = 1e-8;
        maxStepNum = 1e2;
        stepNum
        stepShrnk = 0.8;
        preG=1;
        preP=0;
        cost
        deltaNormAlpha
        PsitPhitz;
        PsitPhit1;
        converged = false;
        warned = false;
    end
    methods
        function obj = ConjGrad(n,alpha)
            if(nargin>0)
                obj.n = n;
            end
            obj.fArray = cell(obj.n,1);
            obj.fVal = zeros(obj.n,1);
            obj.hArray = cell(obj.n,1);
            obj.coef = zeros(obj.n,1);
            obj.alpha = alpha;
        end

        function [f,g,h] = func(obj,aaa)
            if(nargout==3)
                g = 0;
                for i = 1:obj.n
                    [obj.fVal(i),gtemp,obj.hArray{i}] = obj.fArray{i}(aaa);
                    g = g+gtemp*obj.coef(i);
                end
                f = obj.fVal'*obj.coef;
                h = @(xxx,opt) hessian(xxx,opt);
            elseif(nargout==2)
                g = 0;
                for i = 1:obj.n
                    [obj.fVal(i),gtemp] = obj.fArray{i}(aaa);
                    g = g+gtemp*obj.coef(i);
                end
                f = obj.fVal'*obj.coef;
            else
                for i = 1:obj.n
                    obj.fVal(i) = obj.fArray{i}(aaa);
                end
                f = obj.fVal'*obj.coef;
            end
            function hh = hessian(x,opt)
                hh=0;
                for j = 1:obj.n
                    hh = hh + obj.hArray{j}(x,opt)*obj.coef(j);
                end
                if(hh<=0)
                    warning(['fAlpha: obj function is non-convex over alpha,'...
                        '%g'],hh);
                    keyboard
                    obj.warned = true;
                end

            end
        end

        function prCG(obj)
            pp=0; obj.converged = false; obj.warned = false; needBreak = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [oldCost,grad,hessian] = obj.func(obj.alpha);

                maxStep=1;
                beta=grad'*(grad-obj.preG)/(obj.preG'*obj.preG);
                deltaAlpha=grad+max(beta,0)*obj.preP;
                obj.deltaNormAlpha=grad'*deltaAlpha;
                s1=obj.deltaNormAlpha/hessian(deltaAlpha,2);
                %atHessianA(deltaAlpha,weight,t1*hphi,Phi,Phit,t3, Psit,opt.muLustig,sqrtSSqrMu);
                obj.preP = deltaAlpha; obj.preG = grad;
                deltaAlpha = deltaAlpha*s1;
                obj.deltaNormAlpha = obj.deltaNormAlpha*s1;

                if(obj.deltaNormAlpha<obj.thresh)
                    obj.converged=true;
                    break;
                end

                % start of line Search
                ppp=0; stepSz=min(1,maxStep);
                while(true)
                    ppp=ppp+1;
                    newX=obj.alpha-stepSz*deltaAlpha;
                    %newX(newX<0)=0; % force it be positive;
                    newCost=obj.func(newX);

                    if(newCost <= oldCost - stepSz/2*obj.deltaNormAlpha)
                        obj.alpha = newX;
                        obj.cost = newCost;
                        break;
                    else
                        if(ppp==9) keyboard
                        end
                        if(ppp>10)
                            obj.cost = oldCost;
                            obj.warned = true;
                            needBreak = true;
                            warning('exit iterations for higher convergence criteria: %g\n',obj.deltaNormAlpha);
                            break;
                        else
                            stepSz=stepSz*obj.stepShrnk;
                        end
                    end
                end
                %end of line search
                if(obj.converged || needBreak) break; end
            end
            obj.stepNum = pp;
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
                    t1=min([1, abs(diff0'*difphi/norm(difphi)), abs(costA/costB)]);
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

