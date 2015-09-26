classdef AT < Methods
    properties
        stepShrnk = 0.5;
        zbar=0;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        admmAbsTol=1e-9;
        admmTol=1e-3;
        cumu=0;
        cumuTol=4;
        nonInc=0;
        innerSearch=0;

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;

        maxInnerItr=100;

        proxmapping
    end
    methods
        function obj = AT(n,alpha,maxAlphaSteps,stepShrnk,pm)
        %   alpha(alpha<0)=0;
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.nonInc=0;
            obj.zbar=alpha;
            obj.proxmapping=pm;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.cumu=0;
            obj.theta=0;
            obj.zbar=alpha;
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.warned = false;
            pp=0; obj.debug='';
            if(obj.restart>0) obj.restart=0; end

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;
                obj.theta=(1+sqrt(1+4*obj.theta^2))/2;

                y=(1-1/obj.theta)*obj.alpha+obj.zbar/obj.theta;

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

                    [newZbar,obj.innerSearch]=obj.proxmapping(zbar-obj.grad/obj.t*obj.theta,...
                        obj.theta*obj.u/obj.t,obj.admmTol*obj.difAlpha,obj.maxInnerItr);
                    newX=(1-1/obj.theta)*obj.alpha+newZbar/obj.theta;

                    newCost=obj.func(newX);
                    LMM=(oldCost+innerProd(obj.grad,newX-y)+sqrNorm(newX-y)*obj.t/2);
                    if(newCost<=LMM)
                        if(obj.p<=obj.preSteps && obj.ppp<18 && goodStep && obj.t>0)
                            obj.t=obj.t*obj.stepShrnk; continue;
                        else
                            break;
                        end
                    else
                        if(obj.ppp<=20 && obj.t>0)
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
                obj.fVal(3) = obj.fArray{3}(newX);
                temp = newCost+obj.u*obj.fVal(3);

                % restart
                if((temp-obj.cost)>0)
                    if(goodMM)
                        if(sum(abs(xbar-obj.alpha))~=0) % if has monmentum term, restart
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
                                obj.t=obj.t/obj.stepShrnk; obj.cumu=0;
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
                obj.zbar=newZbar;

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
            obj.theta=0; obj.zbar=obj.alpha;
            recoverT=obj.stepSizeInit('hessian');
            obj.t=min([obj.t;max(recoverT)]);
        end
    end
end

