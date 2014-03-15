classdef ADMM_L1 < Methods
    properties
        stepShrnk = 0.5;
        maxItr=1e2;
        preAlpha=0;
        preG=[];
        preY=[];
        s
        Psi_s
        pa
        rho = 1;
        y1=0;
        main
        absTol=1e-4;
    end
    methods
        function obj = ADMM_L1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.coef(1) = 1;
            obj.maxStepNum = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi;
            obj.Psit = Psit;
            obj.s = obj.Psit(alpha);
            obj.Psi_s = obj.Psi(obj.s);
            fprintf('use ADMM_L1 method\n');
            obj.main = @obj.main_0;
            figure(123); figure(386);
        end
        function main_0(obj)
            obj.p = obj.p+1; obj.warned = false;
            %if(obj.rho<10) obj.rho = obj.rho*1.1; end
            %if(obj.absTol>1e-10) obj.absTol=obj.absTol*0.5; end
            subProb = FISTA_NN(2,obj.alpha);
            subProb.absTol = 1e-14;
            subProb.fArray{1} = obj.fArray{1};
            subProb.fArray{2} = @(aaa) Utils.augLag(aaa,obj.Psi_s-obj.y1);
            subProb.coef = [1; obj.rho];
            obj.alpha = subProb.main();

            obj.s = Utils.softThresh(...
                obj.Psit(obj.alpha+obj.y1),...
                obj.u/(obj.rho));
            obj.Psi_s = obj.Psi(obj.s);

            obj.y1 = obj.y1 - (obj.Psi_s-obj.alpha);
            set(0,'CurrentFigure',386);
            semilogy(obj.p,norm(obj.Psi_s-obj.alpha),'.'); hold on;

            obj.func(obj.alpha);
            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(1:obj.n)'*obj.coef(1:obj.n)+obj.u*obj.fVal(obj.n+1);
        end
        function main_1(obj)
            obj.p = obj.p+1; obj.warned = false;

            %y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
            y=obj.alpha;
            obj.preAlpha = obj.alpha;

            if(isempty(obj.preG))
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
                newX = y + (extra)/(obj.t+obj.rho);
                newCost=obj.func(newX);
                if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2)
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end
            obj.alpha = newX;

            obj.s = Utils.softThresh(...
                obj.Psit(obj.alpha+obj.y1),...
                obj.u/(obj.rho));
            obj.Psi_s = obj.Psi(obj.s);

            obj.y1 = obj.y1 - (obj.Psi_s-obj.alpha);

            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(:)'*obj.coef(:);
        end
    end
end

