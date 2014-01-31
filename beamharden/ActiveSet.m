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
        maxStep = 1e3; % default max # of iterations
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

        function init(obj,Ie);
            if(obj.B(end,:)*Ie<obj.b(end)) 
                Ie=obj.b(end)/(obj.B(end,:)*Ie)*Ie;
            end
            obj.Ie=Ie;
            obj.Q = (obj.B*Ie-obj.b<=obj.epsilon);
            obj.Z = null(obj.B(obj.Q,:),'r');
        end

        function out = main(obj)
            pp=0;
            while(pp<maxPP)
                pp=pp+1;
                [cost,grad,hessian] = obj.func(obj.Ie);
                ppp=0;
                while(1)
                    ppp=ppp+1;
                    zhz=Z'*hessian*Z; temp=min(eig(zhz));
                    if(temp<eps)
                        if(-temp>minZHZ) minZHZ=-temp; end
                        zhz=zhz+minZHZ*eye(size(zhz));
                    end

                    deltaIe=Z*inv(zhz)*(Z'*grad);
                    nextGrad=grad-hessian*deltaIe;
                    temp=inv(B(Q,:)*B(Q,:)')*B(Q,:)*nextGrad;
                    lambda=ones(size(Q)); lambda(Q)=temp;

                    %lambda=0;

                    if(any(lambda<0))
                        temp=find(lambda<0);
                        %temp=find(lambda==min(lambda));
                        [~,temp1]=sort(abs(temp-1/2-E/2),'descend');
                        Q(temp(temp1(end)))=false; fprintf('%s\n', char(Q(:)'+'0') );
                        Z = null(B(Q,:),'r');
                        asIdx=asIdx+1;
                        ASactive.itr(asIdx)=p;
                        ASactive.Ie{asIdx}=Ie;
                    else
                        % determine the maximum possible step size
                        constrainMargin = B*Ie-b;
                        if(any(constrainMargin<-eps))
                            display(Ie);
                            %error('current Ie violate B*I>=b constraints');
                        end
                        temp = B*deltaIe;
                        temp1 = inf*ones(size(temp));
                        temp1(temp>eps) = constrainMargin(temp>eps)./temp(temp>eps);
                        maxStep = min( temp1 );
                        temp = find((temp>eps) & (temp1==maxStep) & (~Q));
                        [~,temp1]=sort(abs(temp-1/2-E/2),'descend');
                        collide = zeros(size(Q))==1; collide(temp(temp1(1))) = true;
                        if(maxStep<eps)
                            % if maxStep ==0 find the one with largest temp
                            % use b1 spline will have better performance.
                            if(any(collide) && ppp>50)
                                Q = (Q | collide);
                                fprintf('%s\n', char(Q(:)'+'0') );
                                Z = null(B(Q,:),'r');
                            else
                                IeReady=true;
                                break;
                            end
                        else
                            break;
                        end
                    end
                end
            end
            out.IeSteps(p)=pp;
        end

        function out = main1(obj,Ie)

            % begin line search
            ppp=0; stepSz=min(1,maxStep);
            while(1)
                ppp=ppp+1;
                newX=Ie-stepSz*deltaIe;
                %newX(newX<0)=0; % force it be positive;

                deltaNormIe=grad'*deltaIe;
                [newCost,zmf]=llI(A,newX);

                if(newCost <= cost - stepSz/2*grad'*deltaIe)
                    break;
                else
                    if(pp>10)
                        stepSz=0; else stepSz=stepSz*stepShrnk;
                    end
                end
            end
            % end of line search
            Ie = newX;

            Ie(Ie<eps)=0;
            if(B(end,:)*Ie<b(end))
                Ie=b(end)/(B(end,:)*Ie)*Ie;
            end
            if(stepSz==maxStep)
                Q = (Q | collide);
                fprintf( '%s\n', char(Q(:)'+'0') );
                Z = null(B(Q,:),'r');
                asIdx=asIdx+1;
                ASactive.itr(asIdx)=p;
                ASactive.Ie{asIdx}=Ie;
            end

            out.llI(p)=newCost;




            %if(deltaNormIe<thresh2) IeReady=true; break; end
            if(deltaNormIe<thresh2) IeReady=true; return; end




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
