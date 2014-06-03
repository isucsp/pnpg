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
        grad;

        debug = true;

        restart=0;   % make this value negative to disable restart
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
            if(obj.restart>=0) obj.restart=0; end
            pp=0;
            while(pp<obj.maxItr)
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp;

                % y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
                % y=obj.alpha;
                obj.preAlpha = obj.alpha;

                % [oldCost,obj.grad] = obj.func(y);
                % obj.t = abs(innerProd(obj.grad-obj.preG, y-obj.preY))/...
                %     sqrNorm(y-obj.preY);
                % obj.preG = obj.grad; obj.preY = y;
                [oldCost,obj.grad] = obj.func(y);
                if(obj.debug) figure; end

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp = obj.ppp+1;
                    newX = y - (obj.grad)/(obj.t);
                    newX = obj.innerADMM_v5(newX,obj.t,obj.u,...
                        max(obj.admmTol*obj.difAlpha,obj.admmAbsTol));
                    newCost=obj.func(newX);
                    if(obj.ppp>20 || newCost<=oldCost+innerProd(obj.grad, newX-y)+sqrNorm(newX-y)*obj.t/2)
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                obj.stepSize=1/obj.t;
                if(obj.ppp==1)
                    obj.cumu=obj.cumu+1;
                    if(obj.cumu>=obj.cumuTol)
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                else obj.cumu=0;
                end
                obj.difAlpha=pNorm(newX-obj.alpha)/pNorm(newX);

                obj.fVal(3) = pNorm(obj.Psit(newX),1);
                temp = newCost+obj.u*obj.fVal(3);

                % restart
                if((isempty(obj.cost) || temp>obj.cost) && obj.restart>=0)
                    obj.theta=0; pp=pp-1;
                    if(obj.restart>0) obj.t=obj.t/obj.stepShrnk; end
                    obj.restart=obj.restart+1;
                    if(obj.restart<10) continue; end
                end
                obj.alpha = newX;
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
                difPsi_s=pNorm(temp-Psi_s)/pNorm(temp);

                temp = alpha;
                alpha = (newX+rho*(Psi_s-y1))/(1+rho);
                alpha(alpha<0)=0;
                difAlpha = pNorm(temp-alpha)/pNorm(temp);

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
            if(obj.debug) absTol=1e-16; end
            %while(pppp<1)
            while(true)
                pppp=pppp+1;
                temp = alpha;
                alpha = (newX+rho*(p-y1));
                alpha = obj.Psi(Utils.softThresh(obj.Psit(alpha),u/t))/(1+rho);
                difAlpha = pNorm(temp-alpha)/pNorm(temp);
                if(isnan(difAlpha) && (pNorm(temp)==0) && (pNorm(temp-alpha)==0))
                            difAlpha=0;
                end

                temp=p; p=alpha+y1; p(p<0)=0;
                difP=pNorm(temp-p)/pNorm(temp);
                if(isnan(difP) && (sqrNorm(temp)==0) && (sqrNorm(temp-p)==0))
                    difP=0;
                end

                y1 = y1 +alpha-p;

                if(obj.debug)
                    da(pppp)=difAlpha;
                    dp(pppp)=difP;
                    dap(pppp)=pNorm(alpha-p);
                    ny(pppp) = pNorm(y1);
                end

                if(difAlpha<absTol && difP<absTol) break; end
                if(obj.debug)
                    if(pppp>4000) break; end
                end
            end
            if(obj.debug)
                semilogy(da,'r'); hold on;
                semilogy(dp,'g'); semilogy(dap,'b'); semilogy(ny,'k');
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
            co=sqrNorm(newX-x)/2*obj.t;
            co=co+obj.u*pNorm(obj.Psit(x));
            if(any(x<0)) co=inf; end
        end
    end
end

