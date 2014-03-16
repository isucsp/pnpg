classdef FISTA_NN < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
        maxItr=1e3;
        absTol=1e-4;
    end
    methods
        function obj = FISTA_NN(n,alpha,maxItr,stepShrnk)
            obj = obj@Methods(n,alpha);
            if(nargin>=3) obj.maxItr=maxItr; end
            if(nargin>=4) obj.stepShrnk=stepShrnk; end
        end
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            for pp=1:obj.maxItr
                y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
                obj.preAlpha = obj.alpha;

                if(isempty(obj.preG))
                    [oldCost,grad,hessian] = obj.func(y);
                    obj.t = hessian(grad,2)/(grad'*grad);
                else
                    [oldCost,grad] = obj.func(y);
                    obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                        ((y-obj.preY)'*(y-obj.preY)));
                end
                obj.preG = grad; obj.preY = y;

                % start of line Search
                obj.ppp=0;
                while(1)
                    obj.ppp=obj.ppp+1;
                    newX = y - grad/obj.t;
                    newX(newX<0)=0;
                    newCost=obj.func(newX);
                    if(newCost<=oldCost+grad'*(newX-y)+obj.t/2*(norm(newX-y)^2))
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                obj.difAlpha=norm(newX-obj.alpha).^2/length(newX);
                obj.alpha = newX;
                set(0,'CurrentFigure',123);
                subplot(2,1,1); semilogy(pp,newCost,'.'); hold on;
                subplot(2,1,2); semilogy(pp,obj.difAlpha,'.'); hold on;
                drawnow;
                if(obj.difAlpha<obj.absTol) break; end
            end
            obj.cost = obj.fVal(1:obj.n)'*obj.coef(1:obj.n);
            out = obj.alpha;
        end
    end
end

