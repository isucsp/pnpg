classdef ConjGrad < handle
    properties
        % for common purpose
        mainFunc = @(o) NCG_PR(o);
        method = 'NCG_PR';
        alpha
        n = 3;
        fArray
        fVal
        coef
        hArray
        difAlpha
        difObj
        thresh = 1e-8;
        maxStepNum = 1e2;
        stepNum
        stepShrnk = 0.8;
        preG=1;
        preP=0;
        converged = false;
        warned = false;
        % for NCG_PR

        % for SpaRSA
        Psi     % Psi and Psit are inv dwt transformations
        Psit
        sigma =0.01;    %
        preAlpha = 0;   
        t=-1;   % suppose be larger than L, largest eigenvalue of Hessian
        u       % coefficient of l1 norm term
        M = 5;  % keep the record of objective from the last M iterations
        cost
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
            obj.cost = zeros(obj.M,1);
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
        function set.method(obj,method)
            switch lower(method)
                case lower('sparsa')
                    obj.mainFunc=@(o) o.SpaRSA();
                    obj.method = 'SpaRSA';
                    fprintf('use SpaRSA method\n');
                case lower('NCG_PR')
                    obj.mainFunc=@(o) o.NCG_PR();
                    obj.method = 'NCG_PR';
                    fprintf('use NCG_PR method\n');
                otherwise
                    error('input method for alpha not found\n');
            end
        end
        function set.M(obj,M)
            obj.M = M;
            obj.cost = zeros(obj.M,1);
        end

        function main(obj)
            obj.mainFunc(obj);
        end

        function SpaRSA(obj)
            pp=0; obj.converged = false; obj.warned = false; needBreak = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [oldCost,grad] = obj.func(obj.alpha);
                si = obj.Psit(obj.alpha);
                dsi = obj.Psit(grad);
                obj.t = (grad-obj.preG)'*(grad-obj.preG)/...
                    ((obj.alpha-obj.preAlpha)'*(obj.alpha-obj.preAlpha));
                obj.cost = [oldCost; obj.cost(1:obj.M-1)];

                % start of line Search
                ppp=0;
                while(~obj.converged && ~needBreak)
                    ppp=ppp+1;
                    wi=si-dsi/obj.t;
                    newSi=ConjGrad.softThresh(wi,obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    newCost=obj.func(newX);
                    obj.difAlpha = (newX-obj.alpha)'*(newX-obj.alpha);

                    % the following is the core of SpaRSA method
                    if((newCost <= max(obj.cost) - obj.sigma*obj.t/2*obj.difAlpha)...
                            || (ppp>10))
                        break;
                    else
                        if(ppp>10)
                            warning('exit iterations for higher convergence criteria: %g\n',obj.difAlpha);
                            if(oldCost>=newCost)
                                obj.alpha = newX;
                            else
                                obj.converged = true;
                            end
                            obj.warned = true;
                            needBreak = true;
                        else
                            obj.t=obj.t/obj.stepShrnk;
                        end
                    end
                end
                obj.preG = grad; obj.preAlpha = obj.alpha;
                obj.alpha = newX; obj.difObj = oldCost-newCost;
            end
        end

        function NCG_PR(obj)
            pp=0; obj.converged = false; obj.warned = false; needBreak = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [oldCost,grad,hessian] = obj.func(obj.alpha);

                maxStep=1;
                beta=grad'*(grad-obj.preG)/(obj.preG'*obj.preG);
                deltaAlpha=grad+max(beta,0)*obj.preP;
                deltaNormAlpha=grad'*deltaAlpha;
                s1=deltaNormAlpha/hessian(deltaAlpha,2);
                %atHessianA(deltaAlpha,weight,t1*hphi,Phi,Phit,t3, Psit,opt.muLustig,sqrtSSqrMu);
                obj.preP = deltaAlpha; obj.preG = grad;
                deltaAlpha = deltaAlpha*s1;
                deltaNormAlpha = deltaNormAlpha*s1;

                if(deltaNormAlpha<obj.thresh)
                    obj.converged=true;
                    break;
                end

                % start of line Search
                ppp=0; stepSz=min(1,maxStep);
                while(~obj.converged && ~needBreak)
                    ppp=ppp+1;
                    newX=obj.alpha-stepSz*deltaAlpha;
                    %newX(newX<0)=0; % force it be positive;
                    newCost=obj.func(newX);

                    if((newCost <= oldCost - stepSz/2*deltaNormAlpha)...
                            || (ppp>10 && newCost < oldCost))
                        needBreak = true;
                    else
                        if(ppp>10)
                            warning('exit iterations for higher convergence criteria: %g\n',deltaNormAlpha);
                            obj.warned = true;
                            needBreak = true;
                        else
                            stepSz=stepSz*obj.stepShrnk;
                        end
                    end
                end
                %end of line search

                obj.difAlpha = (newX-obj.alpha)'*(newX-obj.alpha);
                obj.difObj = newCost - oldCost;
                obj.alpha = newX;
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
    methods(Static)
        function y = softThresh(x,thresh)
            idx = abs(x)<=thresh;
            y(x>0) = x(x>0)-thresh;
            y(x<0) = x(x<0)+thresh;
            y(idx) = 0;
        end
    end
end

