classdef PNPG < Methods
    properties

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
    end
end

