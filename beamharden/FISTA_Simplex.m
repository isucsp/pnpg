classdef FISTA_Simplex < handle
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        admmTol=1e-8;
        cumu=0;
        cumuTol=4;
        B
        b
    end
    methods
        function obj = FISTA_Simplex(B,b,alpha,maxItr,stepShrnk)
            obj = obj@Methods(1,alpha);
            obj.B = B; obj.b = b;
            if(nargin>=3) obj.maxItr=maxItr; end
            if(nargin>=4) obj.stepShrnk=stepShrnk; end
        end
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            for pp=1:obj.maxItr
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp;
                obj.preAlpha = obj.alpha;
                %if(isempty(obj.preG))
                %    [oldCost,grad,hessian] = obj.func(y);
                %    obj.t = hessian(grad,2)/(grad'*grad);
                %else
                %    [oldCost,grad] = obj.func(y);
                %    obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                %        ((y-obj.preY)'*(y-obj.preY)));
                %end
                %obj.preG = grad; obj.preY = y;
                if(obj.t==-1)
                    [oldCost,grad,hessian] = obj.func(y);
                    obj.t = hessian(grad,2)/(grad'*grad);
                else
                    [oldCost,grad] = obj.func(y);
                end

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp = obj.ppp+1;
                    newX = y - grad/obj.t;
                    newX(newX<0)=0;
                    newCost=obj.func(newX);
                    if(newCost<=oldCost+grad'*(newX-y)+obj.t/2*(norm(newX-y)^2))
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                if(obj.ppp==1)
                    obj.cumu=obj.cumu+1;
                    if(obj.cumu>=obj.cumuTol)
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                else obj.cumu=0;
                end
                obj.difAlpha=norm(newX-obj.alpha)/norm(newX);
                obj.alpha = newX;
                %if(newCost>obj.cost)
                %    obj.theta=0; obj.preAlpha=obj.alpha;
                %end
                obj.cost = newCost;
                %set(0,'CurrentFigure',123);
                %subplot(2,1,1); semilogy(obj.p,newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                %drawnow;
                if(obj.difAlpha<obj.thresh) break; end
            end
            out = obj.alpha;
        end
        function Ie=adjust(obj,Ie)
            Ie(Ie<0)=0;
            d=obj.b(end)-obj.B(end,:)*Ie;
            while(d>0)
                S = Ie>0 & obj.B(end,:)'<0;
                step = min( min(-Ie(S)./obj.B(end,S)'), ...
                    d/(norm(obj.B(end,S))^2));
                Ie(S) = Ie(S) + step*obj.B(end,S)';
                d=obj.b(end)-obj.B(end,:)*Ie;
            end
        end
    end
end

