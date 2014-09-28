classdef Cell < handle
    methods(Static)
        function summary = getField(ele,field)
            summary=ones(size(ele))*inf;
            for i=1:length(ele(:))
                if(~isempty(ele{i}))
                    res{i}=getfield(ele{i},field);
                    summary(i)=res{i}(end);
                end
            end
        end
    end
end

