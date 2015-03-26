classdef PGs < Methods
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
        adaptiveStep=true;
    end
    methods
        function obj = PGs(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.nonInc=0;
        end
        % solves l(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside IST for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            pp=0; obj.debug='';
            while(pp<obj.maxItr)
                pp=pp+1;
                xbar=obj.alpha;

                [oldCost,obj.grad] = obj.func(xbar);
                si = obj.Psit(xbar); dsi = obj.Psit(obj.grad);

                % start of line Search
                obj.ppp=0; goodStep=true; temp=0; goodMM=true;
                while(true)
                    if(goodStep && temp<obj.adaptiveStep && obj.cumu>=obj.cumuTol)
                        % adaptively increase the step size
                        temp=temp+1;
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                    obj.ppp = obj.ppp+1;
                    wi=si-dsi/obj.t;
                    newSi=Utils.softThresh(wi,obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    obj.newCost=obj.func(newX);
                    LMM=(oldCost+innerProd(obj.grad,newX-xbar)+sqrNorm(newX-xbar)*obj.t/2);
                    if(obj.newCost<=LMM)
                        if(obj.p<=obj.preSteps && obj.ppp<20 && goodStep)
                            obj.t=obj.t*obj.stepShrnk; continue;
                        else
                            break;
                        end
                    else
                        if(obj.ppp<=20)
                            obj.t=obj.t/obj.stepShrnk; goodStep=false; 
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
                if(temp>obj.cost)
                    if(goodMM)
                        obj.debug=[obj.debug 'forceConverge'];
                        newX=obj.alpha;
                        temp=obj.cost;
                    else
                        pp=pp-1;
                        obj.debug=[obj.debug 'falseMonotone'];
                        continue;
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
    end
end

