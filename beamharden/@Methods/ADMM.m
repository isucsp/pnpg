function ADMM(obj)
    pp=0; obj.converged = false; obj.warned = false; needBreak = false;
    while(pp<obj.maxStepNum)
        obj.p = obj.p+1;
        pp=pp+1;

        if(obj.p==1)
            % initialization
            Psi_s = Psi(obj.s);

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
            newX = obj.alpha - grad/t;
            newCost=obj.func(newX);

            if((newCost<=oldCost+obj.t/2*gradSN) || (ppp>10))
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

        obj.alpha = (obj.t*newX+ obj.rho*obj.Psi_s - obj.y1)/(obj.t+obj.rho);
        obj.pa = obj.Psi_s - obj.y2/obj.rho; obj.pa(obj.pa<0) = 0;
        obj.s = soft(...
            Psit(obj.alpha+obj.pa+(obj.y1+obj.y2)/obj.rho)/2,...
            obj.u/obj.t);
        obj.Psi_s = Psi(obj.s);

        obj.y1 = obj.y1 - obj.rho*(obj.Psi_s-obj.alpha);
        obj.y2 = obj.y2 - obj.rho*(obj.Psi_s-obj.pa);

        obj.fVal(2) = sum(abs(obj.s));

        obj.preG = grad; obj.preAlpha = obj.alpha;
        obj.cost = newCost;
    end
end

