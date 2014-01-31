classdef ActiveSet < handle
    % Class help goes here
    properties
        func
        B
        b
        epsilon = 1e-14;
        Q
        Z
        Name % Property help goes here
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
    IeStep = ActiveSet(@(AAA,III) gaussLI(Imea,AAA,III),B,b,Ie);
        function obj = ActiveSet(func, B, b, Ie, epsilon)
            % Method help here
            if(nargin>0)
                obj.func = func; obj.B = B; obj.b = b;
                if(nargin>4) obj.epsilon = epsilon; end
                if(B(end,:)*Ie<b(end)) Ie=b(end)/(B(end,:)*Ie)*Ie; end
                Q = (B*Ie-b<=epsilon);
                Z = null(B(Q,:),'r');
            end
        end

        function out = main()
            ppp=0;
            while(1)
                ppp=ppp+1;
                zhz=Z'*h0*Z; temp=min(eig(zhz));
                if(temp<eps)
                    if(-temp>minZHZ) 
                        minZHZ=-temp;
                    end
                    zhz=zhz+minZHZ*eye(size(zhz));
                end

                deltaIe=zhz\(Z'*diff0); deltaIe=Z*deltaIe;

                estOptG=diff0-h0*deltaIe;
                temp=B(Q,:)*estOptG;
                temp=(B(Q,:)*B(Q,:)')\temp;
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

            % begin line search
            ppp=0; stepSz=min(1,maxStep);
            while(1)
                ppp=ppp+1;
                newX=Ie-stepSz*deltaIe;
                %newX(newX<0)=0; % force it be positive;

                [newCost,zmf]=llI(A,newX);

                if(newCost <= costA - stepSz/2*diff0'*deltaIe)
                    break;
                else
                    if(pp>10)
                        stepSz=0; else stepSz=stepSz*stepShrnk;
                    end
                end
            end
            % end of line search
            deltaNormIe=diff0'*deltaIe;
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
            %if(out.stepSz<1e-10) break; end
            if(deltaNormIe<thresh2)
                IeReady=true; break; end


            end

            function result = backgroundCheck(obj)
                result = queryGovDB(obj.Name,obj.SSNumber);
                if result == false
                    notify(obj,'BackgroundAlert');
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
