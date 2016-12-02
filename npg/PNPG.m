classdef PNPG < Methods
    properties
        stepShrnk = 0.5;
        stepIncre = 0.9;
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

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;

        maxInnerItr=100;
        maxPossibleInnerItr=1e3;

        proximal

        % The following two controls the growth of theta
        gamma=2;
        a=1/4;
    end
    methods
        function obj = PNPG(n,alpha,maxAlphaSteps,stepShrnk,pm)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            if(exist('stepShrnk','var') && ~isempty(stepShrnk))
                obj.stepShrnk = stepShrnk;
            end
            obj.nonInc=0;
            obj.proximal=pm;
            obj.setAlpha(alpha);
            obj.hasLog=true;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.cumu=0;
            obj.theta=0;
            obj.preAlpha=alpha;
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        function out = main(obj)
        end
        function reset(obj)
            obj.theta=1;
            recoverT=obj.stepSizeInit('hessian');
            obj.t=min([obj.t;max(recoverT)]);
        end
    end
end

