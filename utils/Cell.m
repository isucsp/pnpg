classdef Cell < handle
    methods(Static)
        function summary = getField(ele,field,field2)
            summary=ones(size(ele))*inf;
            for i=1:length(ele(:))
                if(~isempty(ele{i}))
                    if(nargin==3)
                        if(isfield(ele{i},field) &&...
                            isfield(ele{i}.(field),field2))
                            res{i}=ele{i}.(field).(field2);
                            summary(i)=res{i}(end);
                        end
                    elseif(nargin==2)
                        if(isfield(ele{i},field))
                            res{i}=ele{i}.(field);
                            summary(i)=res{i}(end);
                        end
                    end
                end
            end
        end
        function summary = getCell(ele,field)
            res=cell(size(ele));
            for i=1:length(ele(:))
                if(~isempty(ele{i}))
                    if(isfield(ele{i},field))
                        res{i}=ele{i}.(field);
                    end
                end
            end
            summary = reshape(res,size(ele));
        end
    end
end

