classdef ActiveSet < handle
    % Find the case when active set wanders
    properties
        func
        B
        b
        epsilon = 1e-15;
        Q
        Z
        Ie
        zmf
        maxItr = 1e2; % default max # of iterations
        thresh = 1e-14; % critera for convergence
        converged = false;
        stepShrnk = 0.8;
        cost = 0;
        stepNum
        course = [];    % record the change of active set
        deltaNormIe
        debugLevel = 0;
        warned = false;
    end 

    methods
        function obj = ActiveSet(B,b,Ie,maxItr,stepShrnk)
            % Method help here
            if(nargin>0)
                obj.B = B; obj.b = b;
                if(nargin>=3)
                    if(nargin>=4)
                        obj.maxItr = maxItr;
                        if(nargin>=5) obj.stepShrnk=stepShrnk; end
                    end
                    obj.init(Ie);
                end
            end
        end

        function init(obj,Ie)
            if(nargin==1) Ie=obj.Ie; end
            obj.Ie = obj.adjust(Ie);
            obj.Q = (obj.B*Ie-obj.b<=eps);
            obj.Z = null(obj.B(obj.Q,:),'r');
        end

        function main(obj)
            obj.converged = false; obj.course = [];
            obj.warned = false; needBreak = false;
            for pp=1:obj.maxItr
                [oldCost,grad,hessian] = obj.func(obj.Ie);
                k = -1*ones(20,1); q = k;
                for ppp=2:20
                    zhz=hessian(obj.Z,3);
                    temp=eig(zhz);
                    maxEig=max(abs(temp)); minEig=min(abs(temp));
                    if((~isempty(temp)) && (min(temp)<0 || (minEig/maxEig)<obj.epsilon))
                        % As this issue is handled gracely, no need to show warning message
                        % if(~obj.warned)
                        %     warning(['\nActiveSet: current Hessian is not',...
                        %         ' positive definite or has bad COND'],0);
                        %     obj.warned = true;
                        % end
                        while(min(temp)<0 || (minEig/maxEig)<obj.epsilon)
                            zhz = zhz + eye(size(zhz))*(2*obj.epsilon*maxEig-minEig)/(1-2*obj.epsilon);
                            temp=eig(zhz);
                            maxEig=max(abs(temp)); minEig=min(abs(temp));
                        end
                    end

                    deltaIe=obj.Z*(zhz\(obj.Z'*grad));
                    nextGrad=grad-hessian(deltaIe,1);
                    temp=(obj.B(obj.Q,:)*obj.B(obj.Q,:)')\(obj.B(obj.Q,:)*nextGrad);
                    lambda=ones(size(obj.Q)); lambda(obj.Q)=temp;

                    if(any(lambda<0))
                        temp=find(lambda<0);
                        [~,temp1]=sort(abs(temp-length(obj.Ie)/2));
                        k(ppp) = temp(temp1(1));
                        %temp=find(lambda==min(lambda));
                        %k(ppp) = temp(randi(length(temp)));
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
                            fprintf('\n'); warning('current Ie violate B*I>=b constraints',0);
                            obj.warned = true;
                        end
                        temp = obj.B*deltaIe;
                        temp1 = inf*ones(size(temp));
                        temp2 = temp>0 & (~obj.Q);
                        temp1(temp2) = constrainMargin(temp2)./temp(temp2);

                        maxStep = min( temp1 );
                        temp = find(temp1==maxStep);
                        [~,temp1]=sort(abs(temp-length(obj.Ie)/2));
                        q(ppp) = temp(temp1(end));
                        collide = zeros(size(obj.Q))==1;
                        collide(q(ppp)) = true;

                        if(maxStep<=0)
                            % if maxStep ==0 find the one with largest temp
                            % use b1 spline will have better performance.
                            if(any(collide) && q(ppp)~=k(ppp-1))
                                break;
                                obj.Q = (obj.Q | collide);
                                obj.Z = null(obj.B(obj.Q,:),'r');
                                obj.course = [obj.course;...
                                    sprintf('%s\n', char(obj.Q(:)'+'0') )];
                            else
                                break;
                            end
                        else
                            break;
                        end
                    end
                end
                if(ppp>=20)
                    fprintf('\n'); warning('Cannot find stable active set, stop at: %s',...
                        sprintf('%s\n', char(obj.Q(:)'+'0') ));
                    obj.warned = true;
                end
                
                obj.deltaNormIe=grad'*deltaIe;

                % begin line search
                ppp=0; stepSz=min(1,maxStep);
                while(~obj.converged && ~needBreak)
                    ppp=ppp+1;
                    newIe=obj.Ie-stepSz*deltaIe;
                    newCost=obj.func(newIe);

                    if((newCost <= oldCost - stepSz/2*obj.deltaNormIe)...
                            || (ppp>10 && newCost < oldCost))
                        obj.Ie = obj.adjust(newIe);
                        obj.cost = newCost;
                        if(stepSz==maxStep)
                            obj.Q = (obj.Q | collide);
                            obj.Z = null(obj.B(obj.Q,:),'r');
                            obj.course = [obj.course;...
                                sprintf('%s\n', char(obj.Q(:)'+'0') )];
                        end      
                        break;
                    else
                        if(ppp>10)
                            fprintf('\n'); warning('exit iterations for higher convergence criteria: %g\n',obj.deltaNormIe);
                            if(oldCost>=newCost)
                                obj.Ie = obj.adjust(newIe);
                                obj.cost = newCost;
                            else
                                fprintf('\n'); warning('ActiveSet: Ie converged\n',0);
                                obj.cost = oldCost;
                                obj.converged = true;
                            end
                            needBreak = true;
                            obj.warned = true;
                        else
                            stepSz=stepSz*obj.stepShrnk;
                        end
                    end
                end
                if(obj.deltaNormIe<obj.thresh)
                    obj.converged=true;
                    break;
                end
                % end of line search
                if(obj.converged || needBreak) break; end
            end
            obj.stepNum = pp;
        end

        function Ie=adjust(obj,Ie)
            Ie(Ie<0)=0;
            if(size(obj.B,1)==2*length(Ie))
                Ie(floor(length(Ie)/2)+1)=max(Ie);
            end
            d=obj.b(end)-obj.B(end,:)*Ie;
            while(d>0)
                S = Ie>0 & obj.B(end,:)'<0;
                step = min( min(-Ie(S)./obj.B(end,S)'), ...
                    d/sqrNorm(obj.B(end,S)));
                Ie(S) = Ie(S) + step*obj.B(end,S)';
                if(size(obj.B,1)==2*length(Ie))
                    Ie(floor(length(Ie)/2)+1)=max(Ie);
                end
                d=obj.b(end)-obj.B(end,:)*Ie;
            end
        end
    end
end
