classdef ADMM_N < Methods
    properties
        s
        Psi_s
        pa
        rho = 0.1;
        y1=0;
    end
    methods
        function obj = ADMM_N(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.coef(1) = 1;
            obj.maxStepNum = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi;
            obj.Psit = Psit;
            obj.s = obj.Psit(alpha);
            obj.Psi_s = obj.Psi(obj.s);
            fprintf('use ADMM_N method\n');
        end
        function main(obj)
            obj.p = obj.p+1; obj.warned = false;

            y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
            %y=obj.alpha;
            obj.preAlpha = obj.alpha;

            if(obj.p==1)
                [oldCost,grad,hessian] = obj.func(y);
                obj.t = hessian(grad,2)/(grad'*grad);
            else
                [oldCost,grad] = obj.func(y);
                obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                    ((y-obj.preY)'*(y-obj.preY)));
            end
            obj.preY = y; obj.preG = grad;
            extra = obj.rho*(obj.Psi_s-obj.y1);

            % start of line Search
            while(1)
                newX = y + (-obj.rho*y-grad+extra)/(obj.t+obj.rho);
                newCost=obj.func(newX);
                if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2)
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end
            obj.alpha = newX;

            obj.s = obj.softThresh(...
                obj.Psit(obj.alpha+obj.y1),...
                obj.u/(obj.rho));
            obj.Psi_s = obj.Psi(obj.s);

            obj.y1 = obj.y1 - (obj.Psi_s-obj.alpha);

            obj.fVal(obj.n+1) = sum(abs(obj.s));

            obj.cost = obj.fVal(:)'*obj.coef(:);
        end
    end
end

