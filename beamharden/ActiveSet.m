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
        Name % Property help goes here
        maxStepNum = 1e3; % default max # of iterations
        thresh = 1e-14; % critera for convergence
        converged = false;
        stepShrnk = 0.8;
    end 

    properties (Dependent)
        JobTitle
    end 

    properties (Transient)
        OfficeNumber
    end 

    properties (SetAccess = protected, GetAccess = private)
        EmpNumber
    end 

    events
        BackgroundAlert
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
            obj.Ie=obj.adjust(Ie);
            obj.Q = (obj.B*Ie-obj.b<=obj.epsilon);
            obj.Z = null(obj.B(obj.Q,:),'r');
        end

        function main(obj)
            pp=0; obj.converged = false;
            while(pp<obj.maxStepNum)
                pp=pp+1;
                [cost,grad,hessian] = obj.func(obj.Ie);
                ppp=0;
                while(1)
                    ppp=ppp+1;
                    zhz=obj.Z'*hessian*obj.Z; temp=min(eig(zhz));
                    if(temp<eps)
                        if(-temp>minZHZ) minZHZ=-temp; end
                        zhz=zhz+minZHZ*eye(size(zhz));
                    end

                    deltaIe=obj.Z*(zhz\(obj.Z'*grad));
                    nextGrad=grad-hessian*deltaIe;
                    temp=(obj.B(obj.Q,:)*obj.B(obj.Q,:)')\obj.B(obj.Q,:)*nextGrad;
                    lambda=ones(size(obj.Q)); lambda(obj.Q)=temp;

                    if(any(lambda<0))
                        temp=find(lambda<0);
                        %temp=find(lambda==min(lambda));
                        [~,temp1]=sort(abs(temp-length(lambda)/2),'descend');
                        obj.Q(temp(temp1(end)))=false; fprintf('%s\n', char(obj.Q(:)'+'0') );
                        obj.Z = null(obj.B(obj.Q,:),'r');
                    else
                        % determine the maximum possible step size
                        constrainMargin = obj.B*obj.Ie-obj.b;
                        if(any(constrainMargin<-eps))
                            display(obj.Ie);
                            error('current Ie violate B*I>=b constraints');
                        end
                        temp = obj.B*deltaIe;
                        temp1 = inf*ones(size(temp));
                        temp1(temp>eps) = constrainMargin(temp>eps)./temp(temp>eps);
                        maxStep = min( temp1 );
                        temp = find((temp>eps) & (temp1==maxStep) & (~obj.Q));
                        [~,temp1]=sort(abs(temp-length(obj.Q)/2),'descend');
                        collide = zeros(size(obj.Q))==1;
                        collide(temp(temp1(1))) = true;
                        if(maxStep<eps)
                            % if maxStep ==0 find the one with largest temp
                            % use b1 spline will have better performance.
                            if(any(collide))
                                obj.Q = (obj.Q | collide);
                                fprintf('%s\n', char(obj.Q(:)'+'0') );
                                obj.Z = null(obj.B(obj.Q,:),'r');
                            else
                                break;
                            end
                        else
                            break;
                        end
                    end
                end

                % begin line search
                ppp=0; stepSz=min(1,maxStep);
                while(1)
                    ppp=ppp+1;
                    newIe=obj.Ie-stepSz*deltaIe;

                    deltaNormIe=grad'*deltaIe;
                    newCost=obj.func(newIe);

                    if(newCost <= cost - stepSz/2*deltaNormIe)
                        break;
                    else
                        if(ppp>10)
                            stepSz=0; else stepSz=stepSz*obj.stepShrnk;
                        end
                    end
                end
                % end of line search
                obj.Ie = obj.adjust(newIe);
                if(stepSz==maxStep)


                    obj.Q = (obj.Q | collide);


                    fprintf( '%s\n', char(obj.Q(:)'+'0') );
                    obj.Z = null(obj.B(obj.Q,:),'r');
                end
                if(deltaNormIe<obj.thresh)
                    obj.converged=true;
                    return;
                end
            end
        end

        function Ie=adjust(obj,Ie)
            Ie(Ie<0)=0;
            d=obj.b-obj.B(end,:)*Ie;
            while(d>0)
                S = Ie>0 & obj.B(end,:)'<0;
                step = min( min(-Ie(S)./obj.B(end,S)), ...
                    d/(norm(obj.B(end,S))^2));
                Ie(S) = Ie(S) + step*obj.B(end,S);
                d=obj.b-obj.B(end,:)*Ie;
            end
        end

        function jobt = get.JobTitle(obj)
            jobt = currentJT(obj.EmpNumber);
        end

        function set.OfficeNumber(obj,setvalue)
            if isInUse(setvalue)
                error('Not available')
            else
                obj.OfficeNumber = setvalue;
            end
        end
    end

    methods (Static)
        function num = getEmpNumber
            num = queryDB('LastEmpNumber') + 1;
        end
    end
end
