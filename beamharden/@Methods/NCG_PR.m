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

        obj.alpha = newX; obj.cost = newCost;
        if(obj.converged || needBreak) break; end
    end
    obj.stepNum = pp;
end
