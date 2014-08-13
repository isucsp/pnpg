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
        main;
        cumu=0;
        cumuTol=4;
        newCost;
        nonInc=0;
        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;
    end
    methods
        function obj = FISTA_L1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            fprintf('use FISTA_L1 method\n');
            obj = obj@Methods(n,alpha);
            obj.maxItr = maxAlphaSteps;
            obj.stepShrnk = stepShrnk;
            obj.Psi = Psi; obj.Psit = Psit;
            obj.preAlpha=alpha;
            obj.prePreAlpha=alpha;
            obj.main=@obj.main_FISTA;
        end
        function out = main_FISTA(obj)
            obj.p = obj.p+1; obj.warned = false;
            if(obj.restart>=0) obj.restart=0; end
            pp=0;
            while(pp<obj.maxItr)
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                obj.theta = temp;

                % y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
                % y=obj.alpha;
                obj.preAlpha = obj.alpha;
                si = obj.Psit(y);

                [oldCost,obj.grad] = obj.func(y);
                dsi = obj.Psit(obj.grad);

                % start of line Search
                obj.ppp=0; temp=true; temp1=0;
                while(true)
                    if(temp && temp1<obj.adaptiveStep && obj.cumu>=obj.cumuTol)
                        % adaptively increase the step size
                        temp1=temp1+1;
                        obj.t=obj.t*obj.stepShrnk;
                    end
                    obj.ppp = obj.ppp+1;
                    wi=si-dsi/obj.t;
                    newSi=Utils.softThresh(wi,obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    obj.newCost=obj.func(newX);
                    if(obj.ppp>20 || obj.newCost<=oldCost+innerProd(obj.grad, newX-y)+sqrNorm(newX-y)*obj.t/2)
                        if(temp && obj.p==1)
                            obj.t=obj.t*obj.stepShrnk;
                            continue;
                        else break;
                        end
                    else obj.t=obj.t/obj.stepShrnk; temp=false;
                    end
                end
                obj.fVal(3) = pNorm(newSi,1);
                temp = obj.newCost+obj.u*obj.fVal(3);

                % restart
                if(obj.restart==0 && (~isempty(obj.cost)) && temp>obj.cost)
                    obj.theta=0; pp=pp-1;
                    % Not necessary, decrease is gauranteed
                    % if(obj.restart>0) obj.t=obj.t/obj.stepShrnk; end
                    obj.restart=1;
                    continue;
                end
                if(temp>obj.cost)
                    obj.nonInc=obj.nonInc+1;
                    if(obj.nonInc>5) newX=obj.alpha; end
                end
                obj.cost = temp;
                obj.stepSize = 1/obj.t;
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

        function out = main_DORE(obj)
            obj.p = obj.p+1; obj.warned = false;
            for pp=1:obj.maxItr
                y=obj.alpha; si = obj.Psit(y);

                if(obj.t==-1)
                    [oldCost,obj.grad,hessian] = obj.func(y);
                    obj.t = hessian(obj.grad,2)/(obj.grad'*obj.grad);
                    if(isnan(obj.t)) obj.t=1; end
                else
                    [oldCost,obj.grad] = obj.func(y);
                end
                dsi = obj.Psit(obj.grad);
                oldCost=oldCost+obj.u*sum(abs(si));

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp=obj.ppp+1;
                    newX=y-obj.grad/obj.t;

                    [oldCost1,grad1,hessian1] = obj.func(newX);
                    temp=hessian1(newX-obj.alpha,2);
                    if(temp==0) c1=0; else c1=-grad1'*(newX-obj.alpha)/temp; end
                    alphaBar=newX+c1*(newX-obj.alpha);

                    [oldCost1,grad1,hessian1] = obj.func(alphaBar);
                    temp=hessian1(alphaBar-obj.preAlpha,2);
                    if(temp==0) c2=0; else c2=-grad1'*(alphaBar-obj.preAlpha)/temp; end
                    newX=alphaBar+c2*(alphaBar-obj.preAlpha);

                    newSi=Utils.softThresh(obj.Psit(newX),obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    obj.newCost=obj.func(newX);
                    obj.newCost=obj.newCost+obj.u*sum(abs(newSi));

                    %if(obj.newCost<=oldCost+obj.grad'*(newX-y)+obj.t/2*(norm(newX-y)^2) || obj.ppp>20)
                    if(obj.newCost<oldCost);
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                obj.newCost=obj.newCost+obj.u*sum(abs(newSi));
                if(obj.ppp==1)
                    obj.cumu=obj.cumu+1;
                    if(obj.cumu>=obj.cumuTol)
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                else obj.cumu=0;
                end
                obj.difAlpha=relativeDif(obj.alpha,newX);
                obj.prePreAlpha=obj.preAlpha; obj.preAlpha=obj.alpha; obj.alpha = newX;

                % decide whether to restart
                obj.fVal(3) = sum(abs(newSi));
                temp = obj.newCost+obj.u*obj.fVal(3);
                %if(temp>obj.cost)
                %    obj.theta=0; obj.preAlpha=obj.alpha;
                %end
                obj.cost = temp;
                set(0,'CurrentFigure',123);
                plot(obj.p,c1,'r.'); hold on;
                plot(obj.p,c2,'g.'); hold on;
                %subplot(2,1,1); semilogy(obj.p,obj.newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                if(obj.difAlpha<=obj.thresh) break; end
            end
            out = obj.alpha;
        end
    end
end
