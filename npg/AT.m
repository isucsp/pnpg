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
        admmTol=1e-2;
        cumu=0;
        cumuTol=4;
        incCumuTol=true;
        nonInc=0;
        innerSearch=0;
        restartEvery=100;

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
            obj.proxmapping=pm;
            obj.setAlpha(alpha);
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
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1; pp=pp+1;
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

                    B=obj.t/obj.preT;
                    newTheta=(1+sqrt(1+4*B*obj.theta^2))/2;
                    xbar=(1-1/newTheta)*obj.alpha+obj.zbar/newTheta;

                    [oldCost,obj.grad] = obj.func(xbar);

                    [newZbar,obj.innerSearch]=obj.proxmapping(...
                        obj.zbar-obj.grad/obj.t*newTheta,...
                        newTheta*obj.u/obj.t,...
                        obj.admmTol*obj.difAlpha,...
                        obj.maxInnerItr,obj.zbar);
                    newX=(1-1/newTheta)*obj.alpha+newZbar/newTheta;

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
                                fprintf('\n AT is having a negative step size, do nothing and return!!');
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

                if((newObj-obj.cost)>0 ||...
                        (mod(obj.p,obj.restartEvery)==0 && obj.restart>=0))
                    if(mod(obj.p,obj.restartEvery)==0 && obj.restart>=0) % if has monmentum term, restart
                        obj.theta=0; obj.zbar=obj.alpha;
                        obj.debug=[obj.debug '_Restart'];
                        global strlen
                        fprintf('\t restart');
                        strlen=0;
                        pp=pp-1; continue;
                    else
                        if(~goodMM)
                            obj.debug=[obj.debug '_Reset'];
                            obj.reset();
                        end
                        if(obj.innerSearch<obj.maxInnerItr && obj.admmTol>1e-6)
                            obj.admmTol=obj.admmTol/10;
                            global strlen
                            fprintf('\n decrease admmTol to %g',obj.admmTol);
                            strlen=0;
                            pp=pp-1; continue;
                        elseif(obj.innerSearch>=obj.maxInnerItr && obj.maxInnerItr<1e3)
                            obj.maxInnerItr=obj.maxInnerItr*10;
                            global strlen
                            fprintf('\n increase maxInnerItr to %g',obj.maxInnerItr);
                            strlen=0;
                            pp=pp-1; continue;
                        end
                        % give up and force it to converge
                        obj.debug=[obj.debug '_ForceConverge'];
                        newObj=obj.cost;  newX=obj.alpha;
                end
                obj.theta = newTheta; obj.zbar=newZbar;
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
            obj.theta=0; obj.zbar=obj.alpha;
            recoverT=obj.stepSizeInit('hessian');
            obj.t=min([obj.t;max(recoverT)]);
        end
    end
end

