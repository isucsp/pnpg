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
            obj.preAlpha=alpha;
            obj.Psi=Psi; obj.Psit=Psit;
        end
        function setAlpha(obj,alpha)
            obj.alpha=alpha;
            obj.cumu=0;
            obj.theta=0;
            obj.preAlpha=alpha;
        end
        % solves L(α) + u*||Ψ'*α||_1
        function out = main(obj)
            obj.warned = false;
            pp=0; obj.debug='';

            while(pp<obj.maxItr)
                obj.p = obj.p+1;
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                xbar=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp; obj.preAlpha = obj.alpha;

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
                            goodMM=false;
                            obj.debug=[obj.debug 'falseMM'];
                            break;
                        end
                    end
                end
                obj.stepSize = 1/obj.t;
                obj.fVal(3) = pNorm(newSi,1);
                temp = newCost+obj.u*obj.fVal(3);

                % restart
                if((temp-obj.cost)>0)
                    if(goodMM)
                        if(pNorm(xbar-obj.alpha,1)~=0) % if has monmentum term, restart
                            if(obj.restart>=0)
                                obj.theta=0;
                                obj.debug=[obj.debug 'restart'];
                                pp=pp-1; continue;
                            end
                        else
                            obj.debug=[obj.debug 'goodMM_but_increasedCost'];
                            global strlen
                            fprintf('\n good MM but increased cost, do nothing\n');
                            strlen=0;
                            obj.cumu=0;
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
end
