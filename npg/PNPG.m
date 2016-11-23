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
            obj.stepShrnk = stepShrnk;
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
            pp=0; obj.debug.clearLog();

            while(pp<obj.maxItr)
                obj.p = obj.p+1; pp=pp+1;
                obj.ppp=0; incStep=false; goodMM=true;
                if(obj.adaptiveStep && obj.cumu>=obj.cumuTol)
                    % adaptively increase the step size
                    obj.t=obj.t*obj.stepIncre;
                    obj.cumu=0;
                    incStep=true;
                end
                % start of line Search
                while(true)
                    obj.ppp = obj.ppp+1;

                    B=obj.t/obj.preT;
                    %newTheta=(1+sqrt(1+4*B*obj.theta^2))/2;
                    if(obj.p==1)
                        newTheta=1;
                    else
                        newTheta=1/obj.gamma+sqrt(obj.a+B*obj.theta^2);
                    end
                    xbar=obj.alpha+(obj.theta -1)/newTheta*(obj.alpha-obj.preAlpha);
                    xbar=obj.prj_C(xbar);

                    [oldCost,obj.grad] = obj.func(xbar);

                    if(obj.proximal.iterative)
                        obj.proximal.thresh = obj.admmTol*obj.difAlpha;
                        obj.proximal.maxItr = obj.maxInnerItr;
                    end

                    newX=obj.proximal.denoise(xbar-obj.grad/obj.t, obj.u/obj.t);

                    newCost=obj.func(newX);
                    %if(obj.p<15) keyboard; end
                    if(Utils.majorizationHolds(newX-xbar,newCost,oldCost,[],obj.grad,obj.t))
                        break;
                    else
                        if(obj.ppp<=20 && obj.t>0)
                            obj.t=obj.t/obj.stepShrnk;
                            if(incStep)
                                if(obj.incCumuTol)
                                    obj.cumuTol=obj.cumuTol+4;
                                end
                                incStep=false;
                            end
                        else  % don't know what to do, mark on debug and break
                            if(obj.t<0)
                                error('\n PNPG is having a negative step size, do nothing and return!!');
                            end
                            goodMM=false;
                            obj.debug.appendLog('_FalseMM');
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = obj.proximal.getPenalty();
                newObj = newCost+obj.u*obj.fVal(3);
                objBar = oldCost+obj.u*obj.proximal.getPenalty(xbar);

                if((newObj-obj.cost)>0)
                    if(goodMM && pNorm(xbar-obj.alpha,1)~=0 && obj.restart>=0) % if has monmentum term, restart
                        obj.theta=1;
                        obj.debug.appendLog('_Restart');
                        obj.debug.printWithoutDel(1,'\t restart');
                        pp=pp-1; continue;
                    elseif((~goodMM) || (objBar<newObj))
                        if(~goodMM)
                            obj.debug.appendLog('_Reset');
                            obj.reset();
                        end
                        if(obj.proximal.iterative)
                            if(obj.proximal.steps<obj.maxInnerItr && obj.admmTol>1e-6)
                                obj.admmTol=obj.admmTol/10;
                                obj.debug.printWithoutDel(1,...
                                    sprintf('\n decrease admmTol to %g',obj.admmTol));
                                pp=pp-1; continue;
                                %% IMPORTANT! if not requir too high accuracy
                                %% use 1e3 for maxInnerItr
                            elseif(obj.proximal.steps>=obj.maxInnerItr && obj.maxInnerItr<obj.maxPossibleInnerItr)
                                obj.maxInnerItr=obj.maxInnerItr*10;
                                obj.debug.printWithoutDel(1,...
                                    sprintf('\n increase maxInnerItr to %g',obj.maxInnerItr));
                                pp=pp-1; continue;
                            end
                        end

                        % give up and force it to converge
                        obj.debug.appendLog('_ForceConverge');
                        newObj=obj.cost;  newX=obj.alpha;
                        obj.innerSearch=0;
                    end
                else
                    if(obj.proximal.iterative)
                        obj.proximal.setInit();
                        obj.innerSearch=obj.proximal.steps;
                    else
                        obj.innerSearch=0;
                    end
                end
                obj.theta = newTheta; obj.preAlpha = obj.alpha;
                obj.cost = newObj;
                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;
                obj.preT=obj.t;

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

