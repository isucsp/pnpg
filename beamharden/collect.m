function res = collect(out,field,t)
    
    for i=1:length(out(:))
        temp = out{i}.(field);
        if(nargin==3 && isa(t,'function_handle'))
            temp = t(temp);
        end
        if(nargin==3 && isa(t,'char') && strcmpi(t,'cell'))
            res{i} = temp;
        else
            if(length(temp(:))==1 && ~isstruct(temp))
                res(i) = temp;
            else
                res{i} = temp;
            end
        end
    end
    res = reshape(res,size(out));
end
