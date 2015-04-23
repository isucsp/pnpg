classdef NPG_ind < handle
    properties
        Ie;
        func;
        debug;
        ppp
        p=0;
        t=-1;
        stepSize;
        cost;
        difIe;
        warned=false;
        zmf;

        stepShrnk = 0.5;
        preSteps=10;
        preIe=0;
        preG=[];
        preY=[];
        thresh=1e-5;
        maxItr=1e3;
        theta = 0;
        admmAbsTol=1e-9;
        admmTol=1e-3;   % abs value should be 1e-8
        cumu=0;
        cumuTol=4;
        nonInc=0;
        innerSearch=0;

        restart=0;   % make this value negative to disable restart
        adaptiveStep=true;

        forcePositive=false;
        maxInnerItr=100;

        proximal
    end
    methods
        function obj = NPG_ind(Ie,nonneg,B,b,maxItr,stepShrnk,thresh)
            B=B(:);
            % subject to B'*x <= b
            obj.Ie=Ie;
            obj.preIe=Ie;
            obj.p=0;
            if(nargin>=5) obj.maxItr=maxItr; end
            if(nargin>=6) obj.stepShrnk=stepShrnk; end
            if(nargin>=7) obj.thresh=thresh; end

            if(nonneg && isempty(B))
                obj.proximal=@(xx) max(xx,0);
            elseif(nonneg && ~isempty(B))
                obj.proximal=@(xx) NPG_ind.simplex(xx,B,b);
            elseif(~nonneg && isempty(B))
                obj.proximal=@(xx) xx;
            else
                obj.proximal=@(xx) NPG_ind.BxLessThanb(xx,B,b);
            end
        end
        function setIe(obj,Ie)
            obj.Ie=Ie;
            obj.preIe=Ie;
            obj.cumu=0;
            obj.theta=0;
        end
        % solves L(α) + I(α>=0) + u*||Ψ'*α||_1
        % method No.4 with ADMM inside FISTA for NNL1
        % the order of 2nd and 3rd terms is determined by the ADMM subroutine
        function out = main(obj)
            obj.p = obj.p+1;
            pp=0; obj.debug='';
            if(obj.restart>0) obj.restart=0; end

            while(pp<obj.maxItr)
                pp=pp+1;
                temp=(1+sqrt(1+4*obj.theta^2))/2;
                xbar=obj.Ie+(obj.theta -1)/temp*(obj.Ie-obj.preIe);
                % if(obj.forcePositive) xbar(xbar<0)=0; end
                obj.theta = temp; obj.preIe = obj.Ie;

                if(obj.t==-1)
                    [oldCost,grad,hessian] = obj.func(xbar);
                    obj.t = hessian(grad,2)/(grad'*grad);
                else
                    [oldCost,grad] = obj.func(xbar);
                end

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

                    [newX] = obj.proximal(xbar-grad/obj.t);
                    newCost=obj.func(newX);
                    LMM=(oldCost+innerProd(grad,newX-xbar)+sqrNorm(newX-xbar)*obj.t/2);
                    if((LMM-newCost)>=0)
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

                % restart
                if((newCost-obj.cost)>0)
                    if(goodMM)
                        if(sum(abs(xbar-obj.Ie))~=0) % if has monmentum term, restart
                            obj.theta=0;
                            obj.restart= 1; % make sure only restart once each iteration
                            obj.debug=[obj.debug 'restart'];
                            pp=pp-1; continue;
                        else
                            obj.debug=[obj.debug 'forceConverge'];
                            newX=obj.Ie;  newCost=obj.cost;
                        end
                    else
                        obj.debug=[obj.debug 'falseMonotone'];
                        pp=pp-1; continue;
                    end
                end
                difCost = (newCost-obj.cost)/max(newCost,1e-12);
                obj.cost = newCost;
                obj.difIe = relativeDif(obj.Ie,newX);
                obj.Ie = newX;

                if(obj.ppp==1 && obj.adaptiveStep)
                    obj.cumu=obj.cumu+1;
                else
                    obj.cumu=0;
                end
%               figure(1); semilogy(pp,obj.cost,'.c');hold on; 
%               figure(2); semilogy(pp,obj.difIe,'.c');hold on; 
%               figure(3); semilogy(pp,obj.t,'.c');hold on; 
%               figure(4); plot(obj.Ie,'*-');
                if(obj.difCost<=obj.thresh) break; end
            end
            out = obj.Ie;
        end
        function reset(obj)
            obj.theta=0; obj.preIe=obj.Ie;
        end
    end
    methods(Static)
        function x=BxLessThanb(x,B,b)
            x=x-B*max(0,B'*x-b)/sqrNorm(B);
        end
        function x=simplex(x,B,b)
            x=max(x,0);
            d=B'*x-b;
            while(d>0)
                S = x>0 & B>0;
                step = min( min(x(S)./B(S)), d/sqrNorm(B(S)) );
                x(S) = x(S) - step*B(S);
                d=B'*x-b;
            end
        end
    end
end

