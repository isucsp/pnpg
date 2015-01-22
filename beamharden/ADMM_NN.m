classdef ADMM_NN < Methods
    properties
        stepShrnk = 0.5;
        maxItr=1e2;
        pa
        rho = 0.1;
        y1=0;
        y2=0;
        r_norm
    end
    methods
        function obj = ADMM_NN(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.coef(1) = 1;
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi;
            obj.Psit = Psit;
            obj.pa = alpha;
            obj.pa(obj.pa<0)=0;
            fprintf('use ADMM_NN method\n');
            figure(123); figure(386);
        end
        function main(obj)
            obj.p = obj.p+1; obj.warned = false;
            for i=1:obj.maxItr
                subProb = NPGs(2,obj.alpha,1e3,obj.stepShrnk,obj.Psi,obj.Psit);
                subProb.u=obj.u;
                subProb.absTol = 1e-14;
                subProb.fArray{1} = obj.fArray{1};
                subProb.fArray{2} = @(aaa) Utils.augLag(aaa,obj.pa-obj.y1);
                subProb.coef = [1; obj.rho];
                obj.alpha = subProb.main();

                obj.pa = obj.alpha + obj.y1; obj.pa(obj.pa<0) = 0;

                obj.y1 = obj.y1 - (obj.pa-obj.alpha);
                set(0,'CurrentFigure',386);
            semilogy(obj.p,norm(obj.pa-obj.alpha),'.'); hold on;
                
                obj.r_norm = norm(obj.alpha-obj.pa);
            end
            obj.func(obj.alpha);
            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(1:obj.n)'*obj.coef(1:obj.n)+obj.u*obj.fVal(obj.n+1);
        end
        function main_1(obj)
            obj.p = obj.p+1; obj.warned = false;

            %y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
            y=obj.alpha;
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
            extra = obj.rho*(obj.Psi_s-obj.y1-y)-grad;

            % start of line Search
            obj.ppp=0;
            while(1)
                obj.ppp = obj.ppp+1;
                newX = y + extra/(obj.t+obj.rho);
                newCost=obj.func(newX);
                if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2)
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end
            obj.alpha = newX;

            if(abs(newCost-oldCost)/oldCost<1e+5)
                obj.pa = obj.Psi_s - obj.y2; obj.pa(obj.pa<0) = 0;
                obj.s = obj.softThresh(...
                    obj.Psit(obj.alpha+obj.pa+obj.y1+obj.y2)/2,...
                    obj.u/(2*obj.rho));
                obj.Psi_s = obj.Psi(obj.s);

                obj.y1 = obj.y1 - (obj.Psi_s-obj.alpha);
                obj.y2 = obj.y2 - (obj.Psi_s-obj.pa);
            end

            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(:)'*obj.coef(:);
        end
        function main_2(obj)
            obj.p = obj.p+1; obj.warned = false;

            %y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
            y=obj.alpha;
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
            extra = obj.rho*(obj.pa-y-obj.y1)-grad;

            % start of line Search
            obj.ppp=0;
            while(1)
                obj.ppp = obj.ppp+1;
                newX = y + (extra)/(obj.t+obj.rho);
                newCost=obj.func(newX);
                if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2)
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end
            obj.alpha = newX;

            if(abs(newCost-oldCost)/oldCost<1e-2)
                obj.s = obj.softThresh(...
                    obj.Psit(obj.pa)-obj.y2,...
                    obj.u/(obj.rho));
                obj.Psi_s = obj.Psi(obj.s);
                obj.pa = (obj.alpha+obj.y1+obj.Psi_s+obj.Psi(obj.y2))/2;
                obj.pa(obj.pa<0) = 0;

                obj.y1 = obj.y1 - (obj.pa-obj.alpha);
                obj.y2 = obj.y2 - (obj.Psit(obj.pa)-obj.s);
            end

            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(:)'*obj.coef(:);
        end
        function main_3(obj)
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
            extra = obj.rho*(obj.pa+obj.y1+obj.Psi_sy-y-y)-grad;

            % start of line Search
            obj.ppp=0;
            while(1)
                obj.ppp = obj.ppp+1;
                newX = y + (extra)/(obj.t+obj.rho*2);
                newCost=obj.func(newX);
                if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2)
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end
            obj.alpha = newX;

            if(abs(newCost-oldCost)/oldCost<1e12)
                obj.s = obj.softThresh(...
                    obj.Psit(obj.alpha)-obj.y2,...
                    obj.u/(obj.rho));
                obj.Psi_sy = obj.Psi(obj.s+obj.y2);
                obj.pa = (obj.alpha-obj.y1);
                obj.pa(obj.pa<0) = 0;

                obj.y1 = obj.y1 - (obj.alpha-obj.pa);
                obj.y2 = obj.y2 - (obj.Psit(obj.alpha)-obj.s);
            end

            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(:)'*obj.coef(:);
        end
    end
end

