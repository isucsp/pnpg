classdef Condat < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        admmAbsTol=1e-9;
        admmTol=1e-2;
        cumu=0;
        cumuTol=4;
        incCumuTol=true;
        nonInc=0;
        innerSearch=0;

        forcePositive=false;

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;

        maxInnerItr=100;

        proxmapping

        nbt=0;
    end
    methods
        function obj = Condat(n,alpha,maxAlphaSteps,stepShrnk,pm)
            %   alpha(alpha<0)=0;
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.nonInc=0;
            obj.preAlpha=alpha;
            obj.proxmapping=pm;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.cumu=0;
            obj.theta=0;
            obj.preAlpha=alpha;
            obj.prox=@(x,t) max(min(x,abs(t)),-abs(t));
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        function out = main(obj)
            obj.warned = false;
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;

                [oldCost,obj.grad] = obj.func(obj.alpha);
                obj.preAlpha = obj.alpha;

                ybar=obj.prox(obj.y    +obj.sigma*obj.Psit(obj.alpha),obj.u);
                xbar=max(obj.alpha-obj.tau*(obj.grad+obj.Psi(2*ybar-obj.y)),0);
                obj.alpha=obj.alpha+obj.rho*(xbar-obj.alpha);
                obj.y = obj.y+obj.rho*(ybar-obj.y);

                newCost=obj.func(newX);
                obj.stepSize = obj.lambda;
                obj.fVal(3) = obj.fArray{3}(newX);
                temp = newCost+obj.u*obj.fVal(3);
                obj.cost = temp;

                obj.difAlpha = relativeDif(obj.preAlpha,newX);
                obj.alpha = newX;

                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
    end
end

