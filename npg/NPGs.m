classdef NPGs < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        cumu=0;
        cumuTol=4;
        incCumuTol=true;
        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;
        % the weight that applied to the l1 norm,
        % make it a vector with the same size as alpha to enable
        weight=1;
        proxmapping
    end
    methods
        function obj = NPGs(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi=Psi; obj.Psit=Psit;
            obj.setAlpha(alpha);
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.cumu=0;
            obj.theta=0;
            obj.preAlpha=alpha;
        end
        % solves L(α) + u*||Ψ'*α||_1
        function out = main(obj)
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1; pp=pp+1;
                newTheta=(1+sqrt(1+4*obj.theta^2))/2;
                xbar=obj.alpha+(obj.theta -1)/newTheta*(obj.alpha-obj.preAlpha);

                [oldCost,obj.grad] = obj.func(xbar);
                si = obj.Psit(xbar); dsi = obj.Psit(obj.grad);

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
                    newSi=Utils.softThresh(si-dsi/obj.t,obj.weight*obj.u/obj.t);
                    newX = obj.Psi(newSi);
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
                                fprintf('\n NPG is having a negative step size, do nothing and return!!');
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
                obj.fVal(3) = pNorm(newSi,1);
                newObj = newCost+obj.u*obj.fVal(3);
                objBar = oldCost+obj.u*pNorm(si,1);

                if((newObj-obj.cost)>0)
                    if(goodMM && pNorm(xbar-obj.alpha,1)~=0 && obj.restart>=0) % if has monmentum term, restart
                        obj.theta=1;
                        obj.debug=[obj.debug '_Restart'];
                        global strlen
                        fprintf('\t restart');
                        strlen=0;
                        pp=pp-1; continue;
                    elseif((~goodMM) || (objBar<newObj))
                        % give up and force it to converge
                        obj.debug=[obj.debug '_ForceConverge'];
                        newObj=obj.cost;  newX=obj.alpha;
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
