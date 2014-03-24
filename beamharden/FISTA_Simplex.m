classdef FISTA_Simplex < handle
    properties
        stepShrnk = 0.5;
        preIe=0;
        preG=[];
        preY=[];
        thresh=1e-9;
        maxItr=1e3;
        theta = 0;
        admmTol=1e-8;
        zmf;
        Ie
        difIe
        cost
        func
        p=0;
        ppp
        warned = false;
        t=-1;
        cumu=0;
        cumuTol=4;
        B
        b
    end
    methods
        function obj = FISTA_Simplex(B,b,Ie,maxItr,stepShrnk)
            obj.Ie=Ie;
            obj.B = B; obj.b = b;
            if(nargin>=4) obj.maxItr=maxItr; end
            if(nargin>=5) obj.stepShrnk=stepShrnk; end
            obj.preIe=obj.Ie;
        end
        function out = main(obj)
            obj.p = obj.p+1; obj.warned = false;
            for pp=1:obj.maxItr
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                y=obj.Ie+(obj.theta -1)/temp*(obj.Ie-obj.preIe);
                obj.theta = temp;
                obj.preIe = obj.Ie;
                if(obj.t==-1)
                    [oldCost,grad,hessian] = obj.func(y);
                    obj.t = hessian(grad,2)/(grad'*grad);
                else
                    [oldCost,grad] = obj.func(y);
                end

                % start of line Search
                obj.ppp=0;
                while(true)
                    obj.ppp = obj.ppp+1;
                    newX = obj.adjust(y - grad/obj.t);
                    newCost=obj.func(newX);
                    if(newCost<=oldCost+grad'*(newX-y)+obj.t/2*(norm(newX-y)^2))
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
                obj.difIe=norm(newX-obj.Ie)/norm(newX);
                obj.Ie = newX;
                %if(newCost>obj.cost)
                %    obj.theta=0; obj.preIe=obj.Ie;
                %end
                obj.cost = newCost;
                set(0,'CurrentFigure',123);
                subplot(2,1,1); semilogy(pp,newCost,'.'); hold on;
                subplot(2,1,2); semilogy(pp,obj.difIe,'.'); hold on;
                drawnow;
                if(obj.difIe<obj.thresh)
                    set(0,'CurrentFigure',123);
                    subplot(2,1,1); hold off;
                    subplot(2,1,2); hold off;
                    break;
                end
            end
            out = obj.Ie;
        end
        function Ie=adjust(obj,Ie)
            Ie(Ie<0)=0;
            d=obj.b(end)-obj.B(end,:)*Ie;
            while(d>0)
                S = Ie>0 & obj.B(end,:)'<0;
                step = min( min(-Ie(S)./obj.B(end,S)'), ...
                    d/(norm(obj.B(end,S))^2));
                Ie(S) = Ie(S) + step*obj.B(end,S)';
                d=obj.b(end)-obj.B(end,:)*Ie;
            end
        end
    end
end

