classdef Cell < handle
    methods(Static)
        function summary = getField(ele,field,field2)
            summary=ones(size(ele))*inf;
            for i=1:length(ele(:))
                if(~isempty(ele{i}))
                    if(nargin==3)
                        res{i}=getfield(getfield(ele{i},field),field2);
                        summary(i)=res{i}(end);
                    elseif(nargin==2)
                        res{i}=getfield(ele{i},field);
                        summary(i)=res{i}(end);
                    end
                end
            end
        end
        function summary = getCell(ele,field)
            for i=1:length(ele(:))
                if(~isempty(ele{i}))
                    res{i}=getfield(ele{i},field);
                end
            end
            summary = reshape(res,size(ele));
        end
    end
end

