classdef Debug < handle
    properties
        level
        strlen
        log
    end
    properties (Access = protected)
    end
    properties (Dependent)
    end
    methods
        function obj = Debug(level)
            obj.level=level;
            obj.strlen=0;
            obj.clearLog();
        end
        function printWithoutDel(obj,l,str)
            if(obj.level>=l)
                fprintf(str);
                obj.strlen=0;
            end
        end
        function print(obj,l,str)
            if(obj.level>=l)
                fprintf(str);
                obj.strlen=obj.strlen+length(str);
            end
        end
        function println(obj,l,str)
            if(obj.level>=l)
                if(exist('str','var'))
                    fprintf(str);
                end
                fprintf('\n');
                obj.strlen=0;
            end
        end
        function clear(obj,l)
            if(obj.level>=l)
                if(obj.strlen>0)
                    fprintf(repmat('\b',1,obj.strlen));
                    obj.strlen=0;
                else
                    fprintf('\n');
                end
            end
        end

        function clearLog(obj)
            obj.log='';
        end
        function appendLog(obj,str)
            obj.log=[obj.log str];
        end
    end
end
