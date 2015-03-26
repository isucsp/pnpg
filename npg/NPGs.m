classdef NPGs < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        prePreAlpha=0;
        cumu=0;
        cumuTol=4;
        newCost;
        nonInc=0;
        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;
    end
    methods
        function obj = NPGs(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.preAlpha=alpha;
            obj.prePreAlpha=alpha;
        end
        % solves L(α) + u*||Ψ'*α||_1
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            pp=0; obj.debug='';
            if(obj.restart>0) obj.restart=0; end
            while(pp<obj.maxItr)
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                xbar=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp; obj.preAlpha = obj.alpha;

                [oldCost,obj.grad] = obj.func(xbar);
                si = obj.Psit(xbar); dsi = obj.Psit(obj.grad);

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
                    newSi=Utils.softThresh(si-dsi/obj.t,obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    obj.newCost=obj.func(newX);
                    LMM=(oldCost+innerProd(obj.grad,newX-xbar)+sqrNorm(newX-xbar)*obj.t/2);
                    if((LMM-obj.newCost)>=0)
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
                obj.fVal(3) = pNorm(newSi,1);
                temp = obj.newCost+obj.u*obj.fVal(3);

                % restart
                if((temp-obj.cost)>0)
                    if(goodMM)
                        if(sum(abs(xbar-obj.alpha))~=0) % if has monmentum term, restart
                            obj.theta=0;
                            obj.restart= 1; % make sure only restart once each iteration
                            obj.debug=[obj.debug 'restart'];
                            pp=pp-1; continue;
                        else
                            obj.debug=[obj.debug 'forceConverge'];
                            newX=obj.alpha; temp=obj.cost;
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
