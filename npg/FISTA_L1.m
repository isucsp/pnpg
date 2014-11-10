classdef FISTA_L1 < Methods
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
        function obj = FISTA_L1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.preAlpha=alpha;
            obj.prePreAlpha=alpha;
        end
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            if(obj.restart>0) obj.restart=0; end
            pp=0;
            while(pp<obj.maxItr)
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp; obj.preAlpha = obj.alpha;

                [oldCost,obj.grad] = obj.func(y);
                si = obj.Psit(y); dsi = obj.Psit(obj.grad);

                % start of line Search
                obj.ppp=0; temp=true; temp1=0;
                while(true)
                    if(temp && temp1<obj.adaptiveStep && obj.cumu>=obj.cumuTol)
                        % adaptively increase the step size
                        temp1=temp1+1;
                        obj.cumu=0;
                        obj.t=obj.t*obj.stepShrnk;
                    end
                    obj.ppp = obj.ppp+1;
                    wi=si-dsi/obj.t;
                    newSi=Utils.softThresh(wi,obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    obj.newCost=obj.func(newX);
                    if(obj.ppp>20 || obj.newCost<=oldCost+innerProd(obj.grad, newX-y)+sqrNorm(newX-y)*obj.t/2)
                        if(obj.ppp<=20 && temp && obj.p==1)
                            obj.t=obj.t*obj.stepShrnk;
                            continue;
                        else break;
                        end
                    else obj.t=obj.t/obj.stepShrnk; temp=false; obj.cumuTol=obj.cumuTol+2;
                    end
                end
                obj.fVal(3) = pNorm(newSi,1);
                temp = obj.newCost+obj.u*obj.fVal(3);
                obj.stepSize = 1/obj.t;

                % restart
                if(temp>obj.cost)
                    switch(obj.restart)
                        case 0
                            obj.theta=0;
                            obj.restart= 1; % make sure only restart once each iteration
                            pp=pp-1;
                            continue;
                        case 1
                            newX=obj.alpha;
                            temp=obj.cost;
                            recoverT=obj.stepSizeInit('bb');
                            obj.t=min([obj.t;recoverT(:)]);
                    end
                end
                obj.cost = temp;
                obj.difAlpha = relativeDif(obj.alpha,newX);
                obj.alpha = newX;

                if(obj.ppp==1 && obj.adaptiveStep) obj.cumu=obj.cumu+1;
                else obj.cumu=0; end
                %set(0,'CurrentFigure',123);
                %subplot(2,1,1); semilogy(obj.p,obj.newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
        function reset(obj)
            obj.theta=0; obj.preAlpha=obj.alpha;
        end
    end
end
