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
        newCost;
        nonInc=0;
        innerSearch

        debug = false;

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;
    end
    methods
        function obj = FISTA_ADMM_NNL1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.preAlpha=alpha;
            obj.nonInc=0;
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
                obj.theta = temp; obj.preAlpha = obj.alpha;

                [oldCost,obj.grad] = obj.func(y);

                % start of line Search
                obj.ppp=0; temp=true; temp1=0;
                while(true)
                    if(temp && temp1<obj.adaptiveStep && obj.cumu>=obj.cumuTol)
                        % adaptively increase the step size
                        temp1=temp1+1;
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                    obj.ppp = obj.ppp+1;
                    newX = y - obj.grad/obj.t;
                    [newX,obj.innerSearch] = obj.innerADMM_v4(obj.Psi,obj.Psit,...
                        newX,obj.u/obj.t,obj.admmTol*obj.difAlpha);
                    obj.newCost=obj.func(newX);
                    if(obj.ppp>20 || obj.newCost<=oldCost+innerProd(obj.grad, newX-y)+sqrNorm(newX-y)*obj.t/2)
                        if(temp && obj.p==1)
                            obj.t=obj.t*obj.stepShrnk;
                            continue;
                        else break;
                        end
                    else obj.t=obj.t/obj.stepShrnk; temp=false;
                    end
                end
                obj.fVal(3) = pNorm(obj.Psit(newX),1);
                temp = obj.newCost+obj.u*obj.fVal(3);

                % restart
                if(obj.restart==0 && (~isempty(obj.cost)) && temp>obj.cost)
                    obj.theta=0; pp=pp-1;
                    obj.restart= 1; % make sure only restart once each iteration
                    continue;
                end
                if(temp>obj.cost)
                    obj.nonInc=obj.nonInc+1;
                    if(obj.nonInc>5) newX=obj.alpha; end
                end
                obj.cost = temp;
                obj.stepSize = 1/obj.t;
                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;

                if(obj.ppp==1 && obj.adaptiveStep) obj.cumu=obj.cumu+1;
                else obj.cumu=0; end
                %set(0,'CurrentFigure',123);
                %subplot(2,1,1); semilogy(obj.p,obj.newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
        function reset(obj)
            obj.theta=0; obj.preAlpha=obj.alpha;
        end
    end
    methods(Static)
        function [alpha,pppp] = innerADMM_v4(Psi,Psit,newX,u,absTol)
            % solve 0.5*||α-a||_2 + I(α>=0) + u*||Ψ'*α||_1
            % a is newX;
            % start an ADMM inside the FISTA
            if(~exist('absTol')) absTol=1e-12; end
            alpha=newX; Psi_s=alpha; y1=0; rho=1; pppp=0;
            while(true)
                pppp=pppp+1;

                s = Utils.softThresh(Psit(alpha+y1),u/rho);
                temp=Psi_s; Psi_s = Psi(s);
                difPsi_s=relativeDif(temp,Psi_s);

                temp = alpha;
                alpha = (newX+rho*(Psi_s-y1))/(1+rho);
                alpha(alpha<0)=0;
                difAlpha = relativeDif(temp,alpha);

                y1 = y1 - (Psi_s-alpha);

                %set(0,'CurrentFigure',123);
                %semilogy(pppp,difAlpha,'r.',pppp,difPsi_s,'g.'); hold on;
                %drawnow;

                if(pppp>1e2) break; end
                if(difAlpha<=absTol && difPsi_s<=absTol) break; end
            end
            % end of the ADMM inside the FISTA
        end
        function p = innerADMM_v5(Psi,Psit,newX,u,absTol)
            % solve 0.5*||α-α_0||_2 + u*||Ψ'*α||_1 + I(α>=0) 
            % α_0 is newX;
            % start an ADMM inside the FISTA
            if(~exist('absTol')) absTol=1e-12; end
            alpha=newX; p=alpha; p(p<0)=0; y1=alpha-p;
            rho=1; pppp=0;
            %if(debug) absTol=1e-16; end
            %while(pppp<1)
            while(true)
                pppp=pppp+1;
                temp = alpha;
                alpha = (newX+rho*(p-y1));
                alpha = Psi(Utils.softThresh(Psit(alpha),u))/(1+rho);
                difAlpha = relativeDif(temp,alpha);

                temp=p; p=alpha+y1; p(p<0)=0;
                difP=relativeDif(temp,p);

                y1 = y1 +alpha-p;

                if(debug)
                    da(pppp)=difAlpha;
                    dp(pppp)=difP;
                    dap(pppp)=pNorm(alpha-p);
                    ny(pppp) = pNorm(y1);
                    co(pppp) = 0.5*sqrNorm(newX-p)+u*pNorm(Psit(p),1);
                end

                if(pppp>1e2) break; end
                if(difAlpha<=absTol && difP<=absTol) break; end
            end
            % end of the ADMM inside the FISTA
        end
    end
end

