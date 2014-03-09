function FISTA(obj)
    pp=0; obj.converged = false; obj.warned = false; needBreak = false;
    while(pp<obj.maxStepNum)
        obj.p = obj.p+1;
        pp=pp+1;

        y = obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
        si = obj.Psit(y);

        if(obj.p==1)
            [oldCost,grad,hessian] = obj.func(y);
            deltaNormAlpha=grad'*grad;
            obj.t = hessian(grad,2)/deltaNormAlpha;
        else
            [oldCost,grad] = obj.func(y);
            obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                ((y-obj.preY)'*(y-obj.preY)));
        end
        %oldCost = oldCost+obj.u*sum(abs(si));
        dsi = obj.Psit(grad);

        % start of line Search
        ppp=0;
        while(~obj.converged && ~needBreak)
            ppp=ppp+1;
            wi=si-dsi/obj.t;
            newSi=Methods.softThresh(wi,obj.u/obj.t);
            newX = obj.Psi(newSi);
            difAlpha = (newX-obj.alpha)'*(newX-obj.alpha);

            newCost=obj.func(newX);
            obj.fVal(3) = sum(abs(newSi));
            %newCost = newCost+obj.u*obj.fVal(3);

            % the following is the core of SpaRSA method
            if((newCost <= oldCost + grad'*(newX-y) + obj.t/2*(norm(newX-y)^2))...
                    || (ppp>10))
                break;
            else
                if(ppp>10)
                    warning('exit iterations for higher convergence criteria: %g\n',difAlpha);
                    obj.warned = true;
                    needBreak = true;
                else
                    obj.t=obj.t/obj.stepShrnk;
                end
            end
        end
        obj.preG = grad; obj.preAlpha = obj.alpha; obj.preY = y;
        obj.alpha = newX;
        obj.cost = obj.fVal(:)'*obj.coef(:);
    end
end

