classdef FISTA_ADMM_NNL1 < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        admmAbsTol=1e-9;
        admmTol=1e-3;   % abs value should be 1e-8
        cumu=0;
        cumuTol=4;
    end
    methods
        function obj = FISTA_ADMM_NNL1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            fprintf('use FISTA_ADMM_NNL1 method\n');
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.preAlpha=alpha;
        end
        % solves l(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            for pp=1:obj.maxItr
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp;

                %y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
                %y=obj.alpha;
                obj.preAlpha = obj.alpha;

                %if(isempty(obj.preG))
                %    [oldCost,grad,hessian] = obj.func(y);
                %    obj.t = hessian(grad,2)/(grad'*grad);
                %else
                %    [oldCost,grad] = obj.func(y);
                %    obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                %        ((y-obj.preY)'*(y-obj.preY)));
                %end
                %obj.preG = grad; obj.preY = y;
                if(obj.t==-1)
                    [oldCost,grad,hessian] = obj.func(y);
                    obj.t = hessian(grad,2)/(grad'*grad);
                    if(isnan(obj.t)) obj.t=1; end
                else
                    [oldCost,grad] = obj.func(y);
                end

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp = obj.ppp+1;
                    newX = y - (grad)/(obj.t);
                    newX = obj.innerADMM_v5(newX,obj.t,obj.u,...
                        max(obj.admmTol*obj.difAlpha,obj.admmAbsTol));
                    newCost=obj.func(newX);
                    if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2 || obj.ppp>20)
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                if(obj.ppp==1)
                    obj.cumu=obj.cumu+1;
                    if(obj.cumu>=obj.cumuTol)
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                else obj.cumu=0;
                end
                obj.difAlpha=norm(newX-obj.alpha)/norm(newX);
                obj.alpha = newX;

                obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
                temp = newCost+obj.u*obj.fVal(obj.n+1);

                %if(temp>obj.cost)
                %    obj.theta=0; obj.preAlpha=obj.alpha;
                %end
                obj.cost = temp;
                %set(0,'CurrentFigure',123);
                %subplot(2,1,1); semilogy(obj.p,newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                if(obj.difAlpha<obj.thresh) break; end
            end
            out = obj.alpha;
        end
        function alpha = innerADMM_v4(obj,newX,t,u,absTol)
            % solve 0.5*t*||α-α_0||_2 + I(α>=0) + u*||Ψ'*α||_1
            % which is equivalent to 0.5*||α-α_0||_2 + I(α>=0) + u/t*||Ψ'*α||_1
            % α_0 is newX;
            % start an ADMM inside the FISTA
            alpha=newX; Psi_s=alpha; y1=0; rho=1; pppp=0;
            while(true)
                pppp=pppp+1;

                s = Utils.softThresh(obj.Psit(alpha+y1),u/(rho*t));
                temp=Psi_s; Psi_s = obj.Psi(s);
                difPsi_s=norm(temp-Psi_s)/norm(temp);

                temp = alpha;
                alpha = (newX+rho*(Psi_s-y1))/(1+rho);
                alpha(alpha<0)=0;
                difAlpha = norm(temp-alpha)/norm(temp);

                y1 = y1 - (Psi_s-alpha);

                %set(0,'CurrentFigure',123);
                %semilogy(pppp,difAlpha,'r.',pppp,difPsi_s,'g.'); hold on;
                %drawnow;

                if(difAlpha<absTol && difPsi_s<absTol) break; end
            end
            % end of the ADMM inside the FISTA
        end
        function p = innerADMM_v5(obj,newX,t,u,absTol)
            % solve 0.5*t*||α-α_0||_2 + u*||Ψ'*α||_1 + I(α>=0) 
            % which is equivalent to 0.5*||α-α_0||_2 + u/t*||Ψ'*α||_1 + I(α>=0) 
            % α_0 is newX;
            % start an ADMM inside the FISTA
            alpha=newX; p=alpha; p(p<0)=0; y1=alpha-p;
            rho=1; pppp=0;
            %while(pppp<1)
            while(true)
                pppp=pppp+1;
                temp = alpha;
                alpha = (newX+rho*(p-y1));
                alpha = obj.Psi(Utils.softThresh(obj.Psit(alpha),u/t))/(1+rho);
                difAlpha = norm(temp-alpha)/norm(temp);

                temp=p; p=alpha+y1; p(p<0)=0;
                difP=norm(temp-p)/norm(temp);

                y1 = y1 +alpha-p;

                %set(0,'CurrentFigure',123);
                %semilogy(pppp,difAlpha,'b.',pppp,difP,'c.'); hold on;
                %drawnow;

                if(difAlpha<absTol && difP<absTol) break; end
            end
            % end of the ADMM inside the FISTA
        end
        function newX = innerProjection2(obj,newX,t,u)
            s=obj.Psit(newX);
            s=Utils.softThresh(s,u/t);
            newX=obj.Psi(s);
            newX(newX<0)=0;
        end
        function co = evaluate(obj,newX,x)
            co=norm(newX-x)^2/2*obj.t;
            co=co+obj.u*sum(abs(obj.Psit(x)));
            if(any(x<0)) co=inf; end
        end
    end
end

