classdef NCG_PR < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=1;
        preP=0;
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        admmTol=1e-8;
        cumu=0;
        cumuTol=4;
    end
    methods
        function obj = NCG_PR(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            fprintf('use NCG_PR method\n');
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.preAlpha = alpha;
        end
        function main(obj)
            pp=0; obj.converged = false; obj.warned = false; needBreak = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [oldCost,obj.grad,hessian] = obj.func(obj.alpha);

                maxStep=1;
                beta=obj.grad'*(obj.grad-obj.preG)/(obj.preG'*obj.preG);
                deltaAlpha=obj.grad+max(beta,0)*obj.preP;
                deltaNormAlpha=obj.grad'*deltaAlpha;
                s1=deltaNormAlpha/hessian(deltaAlpha,2);
                %atHessianA(deltaAlpha,weight,t1*hphi,Phi,Phit,t3, Psit,opt.muLustig,sqrtSSqrMu);
                obj.preP = deltaAlpha; obj.preG = obj.grad;
                deltaAlpha = deltaAlpha*s1;
                deltaNormAlpha = deltaNormAlpha*s1;

                if(deltaNormAlpha<obj.thresh)
                    obj.converged=true;
                    break;
                end

                % start of line Search
                obj.ppp=0; stepSz=min(1,maxStep);
                while(~obj.converged && ~needBreak)
                    obj.ppp=obj.ppp+1;
                    newX=obj.alpha-stepSz*deltaAlpha;
                    %newX(newX<0)=0; % force it be positive;
                    newCost=obj.func(newX);

                    if((newCost <= oldCost - stepSz/2*deltaNormAlpha)...
                            || (obj.ppp>10 && newCost < oldCost))
                        needBreak = true;
                    else
                        if(obj.ppp>10)
                            warning('exit iterations for higher convergence criteria: %g\n',deltaNormAlpha);
                            obj.warned = true;
                            needBreak = true;
                        else
                            stepSz=stepSz*obj.stepShrnk;
                        end
                    end
                end
                %end of line search

                obj.stepSize=stepSz;
                obj.alpha = newX; obj.cost = newCost;
                if(obj.converged || needBreak) break; end
            end
            obj.stepNum = pp;
        end
    end
end
