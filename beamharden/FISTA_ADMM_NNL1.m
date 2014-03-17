classdef FISTA_ADMM_NNL1 < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
    end
    methods
        function obj = FISTA_ADMM_NNL1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.coef(1) = 1;
            obj.maxStepNum = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi;
            obj.Psit = Psit;
            fprintf('use FISTA_ADMM_NNL1 method\n');
        end
        % solves l(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;

            y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
            %y=obj.alpha;
            obj.preAlpha = obj.alpha;

            if(isempty(obj.preG))
                [oldCost,grad,hessian] = obj.func(y);
                obj.t = hessian(grad,2)/(grad'*grad);
            else
                [oldCost,grad] = obj.func(y);
                obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                    ((y-obj.preY)'*(y-obj.preY)));
            end
            obj.preY = y; obj.preG = grad;

            % start of line Search
            obj.ppp=0;
            %while(obj.ppp<1)
            while(true)
                obj.ppp = obj.ppp+1;
                newX = y - (grad)/(obj.t);
                newX = obj.innerADMM_v5(newX,obj.t,obj.u,1e-13);

                %newX(newX<0)=0;
                %temp(:,1) = obj.innerADMM_v4(newX,obj.t,obj.u,1e-13);
                %temp(:,2) = obj.innerADMM_v5(newX,obj.t,obj.u,1e-13);
                %temp(:,3) = obj.innerProjection(newX,obj.t,obj.u);
                %temp(:,4) = obj.innerProjection2(newX,obj.t,obj.u);
                %for(i=1:4)
                %    co(i) = obj.evaluate(newX,temp(:,i));
                %end
                %co
                %keyboard                
                %idx=find(co==min(co));
                %newX=temp(:,idx(1));

                newCost=obj.func(newX);
                if(newCost<=oldCost+grad'*(newX-y)+norm(newX-y)^2*obj.t/2)
                    break;
                else obj.t=obj.t/obj.stepShrnk;
                end
            end

            obj.alpha = newX;
            obj.fVal(obj.n+1) = sum(abs(obj.Psit(obj.alpha)));
            obj.cost = obj.fVal(1:obj.n)'*obj.coef(1:obj.n)+obj.u*obj.fVal(obj.n+1);
            out = obj.alpha;
        end
        function alpha = innerADMM_v4(obj,newX,t,u,absTol)
            % solve 0.5*t*||α-α_0||_2 + I(α>=0) + u*||Ψ'*α||_1
            % α_0 is newX;
            % start an ADMM inside the FISTA
            alpha=newX; Psi_s=alpha; y1=0; rho=1; pppp=0;
            while(true)
                pppp=pppp+1;

                s = Utils.softThresh(obj.Psit(alpha+y1),u/(rho));
                temp=Psi_s; Psi_s = obj.Psi(s);
                difPsi_s=norm(temp-Psi_s).^2/length(temp);

                temp = alpha;
                alpha = (t*newX+rho*(Psi_s-y1))/(t+rho);
                alpha(alpha<0)=0;
                difAlpha = norm(temp-alpha).^2/length(temp);

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
            % α_0 is newX;
            % start an ADMM inside the FISTA
            alpha=newX; p=alpha; p(p<0)=0; y1=alpha-p;
            rho=1; pppp=0;
            %while(pppp<1)
            while(true)
                pppp=pppp+1;
                temp = alpha;
                alpha = (t*newX+rho*(p-y1));
                alpha = obj.Psi(Utils.softThresh(obj.Psit(alpha),u))/(t+rho);
                difAlpha = norm(temp-alpha).^2/length(temp);

                temp=p; p=alpha+y1; p(p<0)=0;
                difP=norm(temp-p).^2/length(temp);

                y1 = y1 +alpha-p;

                %set(0,'CurrentFigure',123);
                %semilogy(pppp,difAlpha,'b.',pppp,difP,'c.'); hold on;
                %drawnow;

                if(difAlpha<absTol && difP<absTol) break; end
            end
            % end of the ADMM inside the FISTA
        end
        function newX = innerProjection(obj,newX,t,u)
            %newX(newX<0)=0;
            s=obj.Psit(newX);
            s(s>0)=s(s>0)-u/t;
            s(s<0)=s(s<0)+u/t;
            newX=obj.Psi(s);
            newX(newX<0)=0;
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

