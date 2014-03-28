classdef FISTA_L1 < Methods
    properties
        stepShrnk = 0.5;
        preAlpha=0;
        preG=[];
        prePreAlpha=0;
        prePreG=[];
        preY=[];
        thresh=1e-4;
        maxItr=1e3;
        theta = 0;
        cumu=0;
        cumuTol=4;
        main;
    end
    methods
        function obj = FISTA_L1(n,alpha,maxAlphaSteps,stepShrnk,Psi,Psit)
            %fprintf('use FISTA_L1 method\n');
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
            for pp=1:obj.maxItr
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.alpha+(obj.theta -1)/temp*(obj.alpha-obj.preAlpha);
                test=(obj.theta -1)/temp;
                obj.theta = temp;
                %y=obj.alpha;
                %y=obj.alpha+(obj.p-1)/(obj.p+2)*(obj.alpha-obj.preAlpha);
                obj.preAlpha = obj.alpha;
                si = obj.Psit(y);

                %if(isempty(obj.preG))
                %    [oldCost,grad,hessian] = obj.func(y);
                %    obj.t = hessian(grad,2)/(grad'*grad);
                %else
                %    [oldCost,grad] = obj.func(y);
                %    obj.t = abs( (grad-obj.preG)'*(y-obj.preY)/...
                %        ((y-obj.preY)'*(y-obj.preY)));
                %end
                %obj.preG = grad; obj.preY = y;
                if(obj.t==-1)
                    [oldCost,grad,hessian] = obj.func(y);
                    obj.t = hessian(grad,2)/(grad'*grad);
                    if(isnan(obj.t)) obj.t=1; end
                else
                    [oldCost,grad] = obj.func(y);
                end
                dsi = obj.Psit(grad);

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp=obj.ppp+1;
                    wi=si-dsi/obj.t;
                    newSi=Utils.softThresh(wi,obj.u/obj.t);
                    newX = obj.Psi(newSi);
                    newCost=obj.func(newX);

                    if(newCost<=oldCost+grad'*(newX-y)+obj.t/2*(norm(newX-y)^2) || obj.ppp>20)
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                if(obj.ppp==1)
                    obj.cumu=obj.cumu+1;
                    if(obj.cumu>=obj.cumuTol)
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                else obj.cumu=0;
                end
                obj.difAlpha=norm(newX-obj.alpha)/norm(newX);
                obj.alpha = newX;

                % decide whether to restart
                obj.fVal(obj.n+1) = sum(abs(newSi));
                temp = newCost+obj.u*obj.fVal(obj.n+1);
                %if(temp>obj.cost)
                %    obj.theta=0; obj.preAlpha=obj.alpha;
                %end
                obj.cost = temp;
                set(0,'CurrentFigure',123);
                plot(obj.p,test,'.'); hold on;
                %subplot(2,1,1); semilogy(obj.p,newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                if(obj.difAlpha<obj.thresh) break; end
            end
            out = obj.alpha;
        end

        function out = main_DORE(obj)
            obj.p = obj.p+1; obj.warned = false;
            for pp=1:obj.maxItr
                y=obj.alpha; si = obj.Psit(y);

                if(obj.t==-1)
                    [oldCost,grad,hessian] = obj.func(y);
                    obj.t = hessian(grad,2)/(grad'*grad);
                    if(isnan(obj.t)) obj.t=1; end
                else
                    [oldCost,grad] = obj.func(y);
                end
                dsi = obj.Psit(grad);
                oldCost=oldCost+obj.u*sum(abs(si));

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp=obj.ppp+1;
                    newX=y-grad/obj.t;

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
                    newCost=obj.func(newX);
                    newCost=newCost+obj.u*sum(abs(newSi));

                    %if(newCost<=oldCost+grad'*(newX-y)+obj.t/2*(norm(newX-y)^2) || obj.ppp>20)
                    if(newCost<oldCost);
                        break;
                    else obj.t=obj.t/obj.stepShrnk;
                    end
                end
                newCost=newCost+obj.u*sum(abs(newSi));
                if(obj.ppp==1)
                    obj.cumu=obj.cumu+1;
                    if(obj.cumu>=obj.cumuTol)
                        obj.t=obj.t*obj.stepShrnk;
                        obj.cumu=0;
                    end
                else obj.cumu=0;
                end
                obj.difAlpha=norm(newX-obj.alpha)/norm(newX);
                obj.prePreAlpha=obj.preAlpha; obj.preAlpha=obj.alpha; obj.alpha = newX;

                % decide whether to restart
                obj.fVal(obj.n+1) = sum(abs(newSi));
                temp = newCost+obj.u*obj.fVal(obj.n+1);
                %if(temp>obj.cost)
                %    obj.theta=0; obj.preAlpha=obj.alpha;
                %end
                obj.cost = temp;
                set(0,'CurrentFigure',123);
                plot(obj.p,c1,'r.'); hold on;
                plot(obj.p,c2,'g.'); hold on;
                %subplot(2,1,1); semilogy(obj.p,newCost,'.'); hold on;
                %subplot(2,1,2); semilogy(obj.p,obj.difAlpha,'.'); hold on;
                if(obj.difAlpha<obj.thresh) break; end
            end
            out = obj.alpha;
        end
    end
end
