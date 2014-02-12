classdef ActiveSet < handle
    % Find the case when active set wanders
    properties
        func
        B
        b
        epsilon = 1e-14;
        Q
        Z
        Ie
        maxStepNum = 1e2; % default max # of iterations
        thresh = 1e-14; % critera for convergence
        converged = false;
        stepShrnk = 0.8;
        cost
        stepNum
        course = [];    % record the change of active set
        deltaNormIe
        debugLevel = 7;
    end 

    methods
        function obj = ActiveSet(func, B, b, Ie, epsilon)
            % Method help here
            if(nargin>0)
                obj.func = func; obj.B = B; obj.b = b;
                if(nargin>3)
                    if(nargin>4) obj.epsilon = epsilon; end
                    obj.init(Ie);
                end
            end
        end

        function init(obj,Ie)
            obj.Ie = obj.adjust(Ie);
            obj.Q = (obj.B*Ie-obj.b<=eps);
            obj.Z = null(obj.B(obj.Q,:),'r');
        end

        function main(obj)
            pp=0; obj.converged = false; obj.course = [];
            obj.warned = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [oldCost,grad,hessian] = obj.func(obj.Ie);
                ppp=1;
                k = -1*ones(20,1); q = k;
                while(ppp<20)
                    ppp=ppp+1;
                    zhz=hessian(obj.Z,2);
                    if(min(eig(zhz))<0)
                        keyboard
                    end

                    deltaIe=obj.Z*(zhz\(obj.Z'*grad));
                    nextGrad=grad-hessian(deltaIe,1);
                    temp=(obj.B(obj.Q,:)*obj.B(obj.Q,:)')\obj.B(obj.Q,:)*nextGrad;
                    lambda=ones(size(obj.Q)); lambda(obj.Q)=temp;

                    if(any(lambda<0))
                        temp=find(lambda<0);
                        %temp=find(lambda==min(lambda));
                        [~,temp1]=sort(abs(temp-length(lambda)/2),'descend');
                        k(ppp) = temp(temp1(end));
                        obj.Q(k(ppp))=false;
                        obj.Z = null(obj.B(obj.Q,:),'r');
                        obj.course = [obj.course;...
                            sprintf('%s\n', char(obj.Q(:)'+'0') )];
                    else
                        % determine the maximum possible step size
                        constrainMargin = obj.B*obj.Ie-obj.b;
                        if(any(constrainMargin<-eps))
                            display(obj.Ie);
                            display(constrainMargin);
                            if(obj.debugLevel)
                                keyboard
                            end
                            warning('current Ie violate B*I>=b constraints');
                        end
                        temp = obj.B*deltaIe;
                        temp1 = inf*ones(size(temp));
                        temp1(temp>eps & (~obj.Q)) = ...
                            constrainMargin(temp>eps & (~obj.Q))./temp(temp>eps & (~obj.Q));
                        maxStep = min( temp1 );
                        temp = find(temp1==maxStep);
                        [~,temp1]=sort(abs(temp-length(obj.Q)/2),'descend');
                        q(ppp) = temp(temp1(1));
                        collide = zeros(size(obj.Q))==1;
                        collide(q(ppp)) = true;
                        if(maxStep<eps)
                            % if maxStep ==0 find the one with largest temp
                            % use b1 spline will have better performance.
                            if(any(collide) && q(ppp)~=k(ppp-1))
                                obj.Q = (obj.Q | collide);
                                obj.Z = null(obj.B(obj.Q,:),'r');
                                obj.course = [obj.course;...
                                    sprintf('%s\n', char(obj.Q(:)'+'0') )];
                            else
                                obj.converged=true;
                                break;
                            end
                        else
                            break;
                        end
                    end
                end
                if(ppp>=20)
                    warning('Cannot find stable active set, stop at: %s',...
                        sprintf('%s\n', char(obj.Q(:)'+'0') ));
                end
                
                % begin line search
                ppp=0; stepSz=min(1,maxStep);
                deltaNormIe=grad'*deltaIe;
                if(deltaNormIe<obj.thresh)
                    obj.converged=true;
                end
                while(~obj.converged)
                    ppp=ppp+1;
                    newIe=obj.Ie-stepSz*deltaIe;
                    newCost=obj.func(newIe);

                    if(newCost <= oldCost - stepSz/2*deltaNormIe)
                        obj.Ie = obj.adjust(newIe);
                        obj.cost = newCost;
                        break;
                    else
                        if(ppp>10)
                            obj.cost = oldCost;
                            obj.converged=true;
                            warning('exit iterations for higher convergence criteria: %g\n',deltaNormIe);
                            obj.warned = true;
                        else stepSz=stepSz*obj.stepShrnk;
                        end
                    end
                end
                % end of line search
                if(stepSz==maxStep)
                    obj.Q = (obj.Q | collide);
                    obj.Z = null(obj.B(obj.Q,:),'r');
                    obj.course = [obj.course;...
                        sprintf('%s\n', char(obj.Q(:)'+'0') )];
                end      
                if(obj.converged)
                    break;
                end
            end
            obj.stepNum = pp;
            obj.deltaNormIe = deltaNormIe;
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
