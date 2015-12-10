classdef PG < Methods
    properties
        stepShrnk = 0.5;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        admmAbsTol=1e-9;
        admmTol=1e-2;   % abs value should be 1e-8
        cumu=0;
        cumuTol=4;
        incCumuTol=true;
        nonInc=0;
        innerSearch=0;

        forcePositive=false;
        adaptiveStep=true;

        maxInnerItr=100;

        proxmapping
    end
    methods
        function obj = PG(n,alpha,maxAlphaSteps,stepShrnk,pm)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.nonInc=0;
            obj.proxmapping=pm;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.cumu=0;
            obj.theta=0;
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside IST for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.warned = false;
            pp=0; obj.debug='';
            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;
                xbar=obj.alpha;

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
                            obj.debug=[obj.debug '_falseMM'];
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = obj.fArray{3}(newX);
                temp = newCost+obj.u*obj.fVal(3);

                if((temp-obj.cost)>0)
                    if(goodMM)
                        if(obj.innerSearch<obj.maxInnerItr)
                            obj.difAlpha=0;
                            obj.debug=[obj.debug '_resetDifAlpha'];
                            pp=pp-1; continue;
                        else
                            obj.debug=[obj.debug '_goodMM.but.increasedCost'];
                            global strlen
                            fprintf('\n good MM but increased cost, do nothing\n');
                            strlen=0;
                        end
                    else
                        obj.debug=[obj.debug '_falseMonotone'];
                        if(obj.innerSearch<obj.maxInnerItr)
                            % otherwise do nothing
                            obj.debug=[obj.debug '_resetDifAlpha1'];
                            obj.difAlpha=0;
                            pp=pp-1; continue;
                        end
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
            recoverT=obj.stepSizeInit('hessian');
            obj.t=min([obj.t;max(recoverT)]);
        end
    end
end

