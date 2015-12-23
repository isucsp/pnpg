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
    end
    methods
        function obj = NPG(n,alpha,maxAlphaSteps,stepShrnk,pm)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.nonInc=0;
            obj.proxmapping=pm;
            obj.setAlpha(alpha);
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
                obj.p = obj.p+1; pp=pp+1;
                newTheta=(1+sqrt(1+4*obj.theta^2))/2;
                xbar=obj.alpha+(obj.theta -1)/newTheta*(obj.alpha-obj.preAlpha);
                if(obj.forcePositive)
                    xbar=max(xbar,0);
                end
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
                        obj.u/obj.t,obj.admmTol*obj.difAlpha,...
                        obj.maxInnerItr,obj.alpha);
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
                            if(obj.t<0)
                                global strlen
                                fprintf('\n NPG is having a negative step size, do nothing and return!!\n');
                                strlen=0;
                                return;
                            end
                            goodMM=false;
                            obj.debug=[obj.debug '_FalseMM'];
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = obj.fArray{3}(newX);
                newObj = newCost+obj.u*obj.fVal(3);
                objBar = oldCost+obj.u*obj.fArray{3}(xbar);

                if((newObj-obj.cost)>0)
                    if(goodMM && pNorm(xbar-obj.alpha,1)~=0 && obj.restart>=0) % if has monmentum term, restart
                        obj.theta=1;
                        obj.debug=[obj.debug '_Restart'];
                        pp=pp-1; continue;
                    elseif((~goodMM) || (objBar<newObj))
                        obj.debug=[obj.debug '_Reset'];
                        if(~goodMM) obj.reset(); end
                        if(obj.innerSearch<obj.maxInnerItr && obj.admmTol>1e-6)
                            obj.admmTol=obj.admmTol/10;
                            global strlen
                            fprintf('\n decrease admmTol to %g\n',obj.admmTol);
                            strlen=0;
                            pp=pp-1; continue;
                        elseif(obj.innerSearch>=obj.maxInnerItr && obj.maxInnerItr<1e3)
                            obj.maxInnerItr=obj.maxInnerItr*10;
                            global strlen
                            fprintf('\n increase maxInnerItr to %g\n',obj.maxInnerItr);
                            strlen=0;
                            pp=pp-1; continue;
                        end
                    end
                end
                obj.theta = newTheta; obj.preAlpha = obj.alpha;
                obj.cost = newObj;
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
            obj.theta=1;
            recoverT=obj.stepSizeInit('hessian');
            obj.t=min([obj.t;max(recoverT)]);
        end
    end
end

