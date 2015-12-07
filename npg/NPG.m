classdef NPG < Methods
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
        function obj = NPG(n,alpha,maxAlphaSteps,stepShrnk,pm)
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
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        function out = main(obj)
            obj.warned = false;
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                xbar=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                if(obj.forcePositive)
                    xbar=max(xbar,0);
                end
                obj.theta = temp; obj.preAlpha = obj.alpha;

                [oldCost,obj.grad] = obj.func(xbar);

                obj.ppp=0; goodStep=true; incStep=false; goodMM=true;
                if(obj.adaptiveStep && obj.cumu>=obj.cumuTol)
                    % adaptively increase the step size
                    obj.t=obj.t*obj.stepShrnk;
                    obj.cumu=0;
                    incStep=true;
                end
                % start of line Search
                while(true)
                    obj.ppp = obj.ppp+1;
                    [newX,obj.innerSearch]=obj.proxmapping(xbar-obj.grad/obj.t,...
                        obj.u/obj.t,obj.admmTol*obj.difAlpha,obj.maxInnerItr);
                    newCost=obj.func(newX);
                    if(Utils.majorizationHolds(newX-xbar,newCost,oldCost,[],obj.grad,obj.t))
                        if(obj.p<=obj.preSteps && obj.ppp<18 && goodStep && obj.t>0)
                            obj.t=obj.t*obj.stepShrnk; continue;
                        else
                            break;
                        end
                    else
                        if(obj.ppp<=20 && obj.t>0)
                            obj.t=obj.t/obj.stepShrnk; goodStep=false; 
                            if(incStep)
                                if(obj.incCumuTol)
                                    obj.cumuTol=obj.cumuTol+4;
                                end
                                incStep=false;
                            end
                        else  % don't know what to do, mark on debug and break
                            goodMM=false;
                            obj.debug=[obj.debug 'falseMM'];
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = obj.fArray{3}(newX);
                temp = newCost+obj.u*obj.fVal(3);

                % restart
                if((temp-obj.cost)>0 && obj.restart>=0)
                    if(goodMM)
                        if(sum(abs(xbar-obj.alpha))~=0) % if has monmentum term, restart
                            obj.theta=0;
                            obj.debug=[obj.debug 'restart'];
                            pp=pp-1; continue;
                        else
                            if(obj.innerSearch<obj.maxInnerItr)
                                obj.difAlpha=0;
                                obj.debug=[obj.debug 'resetDifAlpha'];
                                pp=pp-1; continue;
                            else
                                obj.debug=[obj.debug 'goodMM_but_increasedCost'];
                                global strlen
                                fprintf('\n good MM but increased cost, do nothing\n');
                                strlen=0;
                                
%                               obj.t=obj.t/obj.stepShrnk; obj.cumu=0;
%                               newX=obj.alpha;  temp=obj.cost;
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
    methods (Access = protected)
        % this method can be redefined in the subclasses for an indicator
        % of a constraints.
        function res=indicate(obj)
            if(any(obj.alpha<0)) res=inf;
            else res=0; end
        end
    end
    methods(Static)
        function [alpha,pppp] = ADMM(Psi,Psit,a,u,absTol,maxItr,isInDebugMode)
            % solve 0.5*||α-a||_2 + I(α≥0) + u*||Ψ'*α||_1
            if((~exist('absTol','var')) || isempty(absTol)) absTol=1e-6; end
            if((~exist('maxItr','var')) || isempty(maxItr)) maxItr=1e3;  end
            if((~exist('isInDebugMode','var')) || isempty(isInDebugMode)) isInDebugMode=false;  end
            % this makes sure the convergence criteria is nontrival
            absTol=min(1e-3,absTol);
            nu=0; rho=1; cnt=0; preS=Psit(a); s=preS;

            pppp=0;
            while(true)
                pppp=pppp+1;
                cnt= cnt + 1;

                alpha = max((a+rho*Psi(s+nu))/(1+rho),0);
                Psit_alpha=Psit(alpha);
                s = Utils.softThresh(Psit_alpha-nu,u/rho);
                nu=nu+s-Psit_alpha;

                difS=pNorm(s-preS); preS=s;
                residual = pNorm(s-Psit_alpha);
                sNorm = pNorm(s);

                if(isInDebugMode)
                    cost(pppp)=0.5*sqrNorm(max(Psi(s),0)-a)+u*pNorm(Psit(max(Psi(s),0)),1);
                    cost1(pppp)=0.5*sqrNorm(alpha-a)+u*pNorm(Psit_alpha,1);
                    cost2(pppp)=0.5*sqrNorm(alpha-a)+u*pNorm(s,1);

                    rhoRec(pppp)=rho;

                    if(pppp>1)
                        difAlpha = pNorm(preAlpha-alpha);
                        difAlphaRec(pppp-1)=difAlpha/sNorm;
                        difSRec(pppp-1)=difS/sNorm;
                        residualRec(pppp-1)=residual/sNorm;
                    end

                    preAlpha=alpha;
                end

                if(pppp>maxItr) break; end
                if(difS<=absTol*sNorm && residual<=absTol*sNorm) break; end
                if(cnt>10) % prevent excessive back and forth adjusting
                    if(difS>10*residual)
                        rho = rho/2 ; nu=nu*2; cnt=0;
                    elseif(difS<residual/10)
                        rho = rho*2 ; nu=nu/2; cnt=0;
                    end
                end
            end 
            alpha= max(Psi(s),0);

            if(isInDebugMode)
                costRef=0.5*sqrNorm(max(a,0)-a)+u*pNorm(Psit(max(a,0)),1);
                figure;
                semilogy(rhoRec,'r'); title('rho');
                figure;
                semilogy(difAlphaRec,'r'); hold on;
                semilogy(difSRec,'g');
                semilogy(residualRec,'b');
                legend('difAlpha','difS','residual');
                title('covnergence criteria');
                figure;
                semilogy(cost1-min([cost,cost1,cost2]),'r-.'); hold on;
                semilogy(cost -min([cost,cost1,cost2]),'g-'); hold on;
                semilogy(cost2 -min([cost,cost1,cost2]),'b-'); hold on;
                title('admm centered obj');
                legend('alpha','s','alpha,s');
                figure;
                semilogy(cost1,'r'); hold on;
                semilogy(cost,'g'); hold on;
                semilogy(ones(size(cost))*costRef,'k');
                title('admm obj');
                legend('alpha','s','ref');
            end
            % end of the ADMM inside the NPG
        end
    end
end

