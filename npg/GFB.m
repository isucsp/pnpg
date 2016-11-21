classdef GFB < Methods
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

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;

        maxInnerItr=100;

        proxmapping

        nbt=0;
        z1
        z2
        r
        lambda
        w
        prox
    end
    methods
        function obj = GFB(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit,L)
            %   alpha(alpha<0)=0;
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.prox=@(x,t) obj.Psi(Utils.softThresh(obj.Psit(x),t));
            obj.setAlpha(alpha);
            obj.r=1.8/L;
            obj.lambda=1;
            obj.w=0.5;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;

            obj.z1=alpha;
            obj.z2=alpha;
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        function out = main(obj)
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;

                [oldCost,obj.grad] = obj.func(obj.alpha);

                obj.z1=obj.z1+obj.lambda*(obj.prox(2*obj.alpha-obj.z1-obj.r*obj.grad,obj.r*obj.u/obj.w)-obj.alpha);
                obj.z2=obj.z2+obj.lambda*(     max(2*obj.alpha-obj.z2-obj.r*obj.grad,0                )-obj.alpha);
                newX  =obj.w*obj.z1+(1-obj.w)*obj.z2;

                newCost=obj.func(newX);
                obj.stepSize = obj.lambda;
                obj.fVal(3) = obj.fArray{3}(newX);
                temp = newCost+obj.u*obj.fVal(3);
                obj.cost = temp;

                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;

                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
    end
end

