classdef TestObj < handle
    % the result shows that reading from obj is 10 times slower than from
    % normal variable, and writing is even slower.
    %
    % The differences between x and opt.x are negligible
    properties
        x
    end
    methods
        function y = objxAdd(obj,n)
            obj.x=0;
            for i=1:n
                obj.x=obj.x+1;
            end
            y=obj.x;
        end
        function x = xAdd(obj,n)
            x=0;
            for i=1:n
                x=x+1;
            end
            obj.x=x;
        end
        function x = optxAdd(obj,n)
            opt.x=0;
            for i=1:n
                opt.x=opt.x+1;
            end
            obj.x=opt.x;
        end
        function y = objxRead(obj,n)
            obj.x=3;
            y=0;
            for i=1:n
                y=y+obj.x+1;
            end
            obj.x=y;
        end
        function y = xRead(obj,n)
            x=3;
            y=0;
            for i=1:n
                y=y+x+1;
            end
            obj.x=y;
        end
        function y = optxRead(obj,n)
            opt.x=3;
            y=0;
            for i=1:n
                y=y+opt.x+1;
            end
            obj.x=y;
        end
    end
end
