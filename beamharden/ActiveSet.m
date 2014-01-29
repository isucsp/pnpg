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
