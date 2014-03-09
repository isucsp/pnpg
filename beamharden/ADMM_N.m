classdef ADMM_N < Methods
    properties
        s
        Psi_s
        pa
        rho = 0.01;
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
            obj.method = 'ADMM_N';
        end
        function main(obj)
            pp=0; obj.converged = false; obj.warned = false; needBreak = false;
            while(pp<obj.maxStepNum)
                obj.p = obj.p+1;
                pp=pp+1;

                if(obj.p==1)
                    % initialization
                    obj.Psi_s = obj.Psi(obj.s);

                    [oldCost,grad,hessian] = obj.func(obj.alpha);
                    obj.t = hessian(grad,2)/(grad'*grad);
                else
                    [oldCost,grad] = obj.func(obj.alpha);
                    obj.t = abs( (grad-obj.preG)'*(obj.alpha-obj.preAlpha)/...
                        ((obj.alpha-obj.preAlpha)'*(obj.alpha-obj.preAlpha)));
                end
                gradSN = norm(grad)^2;

                % start of line Search
                ppp=0;
                while(~obj.converged && ~needBreak)
                    ppp=ppp+1;
                    newX = obj.alpha - grad/obj.t;
                    newCost=obj.func(newX);

                    if((newCost<=oldCost-gradSN/2/obj.t) || (ppp>10))
                        break;
                    else
                        if(ppp>10)
                            warning('exit iterations for higher convergence criteria: %g\n',difAlpha);
                            obj.warned = true; needBreak = true;
                        else
                            obj.t=obj.t/obj.stepShrnk;
                        end
                    end
                end

                obj.preAlpha = obj.alpha;
                obj.alpha = (obj.t*newX+ obj.rho*obj.Psi_s - obj.y1)/(obj.t+obj.rho);
                obj.s = ADMM.softThresh(...
                    obj.Psit(obj.alpha+obj.y1/obj.rho),...
                    obj.u/(obj.rho));
                obj.Psi_s = obj.Psi(obj.s);

                obj.y1 = obj.y1 - obj.rho*(obj.Psi_s-obj.alpha);

                obj.fVal(obj.n+1) = sum(abs(obj.s));

                obj.preG = grad;
                obj.cost = obj.fVal(:)'*obj.coef(:);
            end
        end
    end
end

