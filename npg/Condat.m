classdef Condat < Methods
    properties
        stepShrnk = 0.5;
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
        prox
        y
        tau
        rho
    end
    methods
        function obj = Condat(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit,L)
            %   alpha(alpha<0)=0;
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.prox=@(x,t) max(min(x,abs(t)),-abs(t));
            obj.setAlpha(alpha);
            obj.sigma=1;
            obj.tau=1/(obj.sigma+L/2) * 0.8;
            obj.rho=(2-L/2/(1/obj.tau-obj.sigma)) * 1;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.y=zeros(size(alpha));
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        function out = main(obj)
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;

                [oldCost,obj.grad] = obj.func(obj.alpha);

                ybar=obj.prox(obj.y    +obj.sigma*obj.Psit(obj.alpha),obj.u);
                xbar=max(obj.alpha-obj.tau*(obj.grad+obj.Psi(2*ybar-obj.y)),0);
                newX=obj.alpha+obj.rho*(xbar-obj.alpha);
                obj.y = obj.y+obj.rho*(ybar-obj.y);

                newCost=obj.func(newX);
                obj.stepSize = obj.rho;
                obj.fVal(3) = obj.fArray{3}(newX);
                obj.cost = newCost+obj.u*obj.fVal(3);

                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;

                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
    end
end

