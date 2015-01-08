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
        innerSearch=0;

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;

        forcePositive=false;
        maxInnerItr=100;
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
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            pp=0; obj.debug='';
            if(obj.restart>0) obj.restart=0; end
            while(pp<obj.maxItr)
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                % if(obj.forcePositive) y(y<0)=0; end
                obj.theta = temp; obj.preAlpha = obj.alpha;

                [oldCost,obj.grad] = obj.func(y);

                % start of line Search
                obj.ppp=0; goodStep=true; incStep=false; goodMM=true;
                while(true)
                    if(obj.adaptiveStep && ~incStep && obj.cumu>=obj.cumuTol)
                        % adaptively increase the step size
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                        incStep=true;
                    end
                    obj.ppp = obj.ppp+1;
                    newX = y - obj.grad/obj.t;
                    [newX,obj.innerSearch] = obj.adaptiveADMM(obj.Psi,obj.Psit,...
                        newX,obj.u/obj.t,obj.admmTol*obj.difAlpha,obj.maxInnerItr);
                    % newX = constrainedl2l1denoise(newX,obj.Psi,obj.Psit,obj.u./obj.t,0,...
                    %      1,1e2,1,obj.admmTol*obj.difAlpha);
                    obj.newCost=obj.func(newX);
                    LMM=(oldCost+innerProd(obj.grad,newX-y)+sqrNorm(newX-y)*obj.t/2);
                    if(obj.newCost<=LMM)
                        if(obj.p<=obj.preSteps && obj.ppp<18 && goodStep)
                            obj.t=obj.t*obj.stepShrnk; continue;
                        else
                            break;
                        end
                    else
                        if(obj.ppp<=20)
                            obj.t=obj.t/obj.stepShrnk; goodStep=false; 
                            if(incStep)
                                obj.cumuTol=obj.cumuTol+4;
                                incStep=false;
                            end
                        else
                            goodMM=false;
                            obj.debug=[obj.debug 'falseMM'];
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = pNorm(obj.Psit(newX),1);
                temp = obj.newCost+obj.u*obj.fVal(3);

                % restart
                if(temp>obj.cost)
                    if(goodMM)
                        if(sum(abs(y-obj.alpha))~=0) % if has monmentum term, restart
                            obj.theta=0;
                            obj.restart= 1; % make sure only restart once each iteration
                            obj.debug=[obj.debug 'restart'];
                            pp=pp-1; continue;
                        else
                            if(obj.innerSearch<obj.maxInnerItr)
                                obj.restart= 2;
                                obj.difAlpha=0;
                                obj.debug=[obj.debug 'resetDifAlpha'];
                                pp=pp-1; continue;
                            else
                                obj.debug=[obj.debug 'forceConverge'];
                                newX=obj.alpha;  temp=obj.cost;
                            end
                        end
                    else
                        obj.debug=[obj.debug 'falseMonotone'];
                        pp=pp-1; continue;
                    end
                end
                obj.cost = temp;
                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;

                if(obj.ppp==1 && obj.adaptiveStep)
                    obj.cumu=obj.cumu+1;
                else
                    obj.cumu=0;
                end
                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
        function reset(obj)
            obj.theta=0; obj.preAlpha=obj.alpha;
            recoverT=obj.stepSizeInit('hessian');
            obj.t=min([obj.t;max(recoverT)]);
        end
    end
    methods(Static)
        function [alpha,pppp] = adaptiveADMM(Psi,Psit,newX,u,absTol,maxItr)
            % solve 0.5*||α-a||_2 + I(α>=0) + u*||Ψ'*α||_1
            % a is newX;
            % start an ADMM inside the FISTA
            if(~exist('absTol','var')) absTol=1e-6; end
            if(~exist('maxItr','var')) maxItr=1e3;  end
            % this makes sure the convergence criteria is nontrival
            absTol=min(1e-3,absTol);
            % if(nargin>5)  figure(123); figure(125);  end
            alpha=newX; Psi_s=alpha; y1=0; rho=1; cnt=0;

            pppp=0;
            while(true)
                pppp=pppp+1;
                cnt= cnt + 1;

                s = Utils.softThresh(Psit(alpha-y1),u/rho);
                temp=Psi_s; Psi_s = Psi(s);
                difPsi_s=pNorm(temp-Psi_s);

                temp = alpha;
                alpha = (newX+rho*(Psi_s+y1))/(1+rho);
                alpha(alpha<0)=0;
                difAlpha = pNorm(temp-alpha);

                y1 = y1 + (Psi_s-alpha);

                residual = pNorm(Psi_s-alpha);
                alphaNorm = pNorm(alpha);

                % if(nargin>5)
                %     set(0,'CurrentFigure',123);
                %     semilogy(pppp,difPsi_s,'r.',pppp,difAlpha,'g.',pppp,residual,'b.'); hold on;
                %     drawnow;
                %     set(0,'CurrentFigure',125);
                %     semilogy(pppp,0.5*sqrNorm(alpha-newX)+u*pNorm(Psit(alpha),1),'g.'); hold on;
                %     drawnow;
                % end

                if(pppp>maxItr) break; end
                if(difAlpha<=absTol*alphaNorm && residual<=absTol*alphaNorm) break; end
                if(cnt>10) % prevent back and forth adjusting
                    if(difAlpha>10*residual)
                        rho = rho/2 ; y1=y1*2; cnt=0;
                    elseif(difAlpha<residual/10)
                        rho = rho*2 ; y1=y1/2; cnt=0;
                    end
                end
            end 
            % end of the ADMM inside the FISTA
        end
        function [alpha,pppp] = innerADMM_v4(Psi,Psit,newX,u,absTol,rho,maxItr)
            % solve 0.5*||α-a||_2 + I(α>=0) + u*||Ψ'*α||_1
            % a is newX;
            % start an ADMM inside the FISTA
            if(~exist('absTol','var')) absTol=1e-12; end
            % this makes sure the convergence criteria is nontrival
            absTol=min(1e-3,absTol);
            alpha=newX; Psi_s=alpha; y1=0;


            if(~exist('rho','var')) rho=1; end
            if(~exist('maxItr','var')) maxItr=1e2; end
            if(nargin>5) figure(123); figure(125); end

            pppp=0;
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

                if(nargin>5)
                    set(0,'CurrentFigure',123);
                    semilogy(pppp,difPsi_s,'r.',pppp,difAlpha,'g.',pppp,relativeDif(Psi_s,alpha),'b.'); hold on;
                    drawnow;
                    set(0,'CurrentFigure',125);
                    semilogy(pppp,0.5*sqrNorm(alpha-newX)+u*pNorm(Psit(alpha),1),'b.'); hold on;
                    drawnow;
                end

                if(pppp>maxItr) break; end
                if(difAlpha<=absTol && difPsi_s<=absTol) break; end
            end
            % end of the ADMM inside the FISTA
        end
        function p = innerADMM_v5(Psi,Psit,newX,u,absTol)
            % solve 0.5*||α-α_0||_2 + u*||Ψ'*α||_1 + I(α>=0) 
            % α_0 is newX;
            % start an ADMM inside the FISTA
            if(~exist('absTol','var')) absTol=1e-12; end
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

