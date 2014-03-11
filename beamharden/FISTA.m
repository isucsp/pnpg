classdef FISTA < Methods
    properties
        s
        Psi_s
    end
    methods
        function obj = FISTA(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.coef(1) = 1;
            obj.maxStepNum = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi;
            obj.Psit = Psit;
            obj.s = obj.Psit(alpha);
            obj.Psi_s = obj.Psi(obj.s);
            fprintf('use FISTA method\n');
        end
        function main(obj)
            obj.p = obj.p+1; obj.warned = false;

            y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
            obj.preAlpha = obj.alpha;
            si = obj.Psit(y);

            if(obj.p==1)
                [oldCost,grad,hessian] = obj.func(y);
                obj.t = hessian(grad,2)/(grad'*grad);
            else
                [oldCost,grad] = obj.func(y);
                obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                    ((y-obj.preY)'*(y-obj.preY)));
            end
            obj.preG = grad; obj.preY = y;
            %oldCost = oldCost+obj.u*sum(abs(si));
            dsi = obj.Psit(grad);

            % start of line Search
            obj.ppp=0;
            while(1)
                obj.ppp=obj.ppp+1;
                wi=si-dsi/obj.t;
                newSi=Methods.softThresh(wi,obj.u/obj.t);
                newX = obj.Psi(newSi);
                newCost=obj.func(newX);

                if(newCost<=oldCost+grad'*(newX-y)+obj.t/2*(norm(newX-y)^2))
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end
            obj.alpha = newX;
            obj.fVal(3) = sum(abs(newSi));
            obj.cost = obj.fVal(:)'*obj.coef(:);
        end
    end
end

