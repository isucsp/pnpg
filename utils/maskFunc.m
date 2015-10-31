
function x = maskFunc(s,maskIdx,n)
% or call this function by maskFunc(s,mask)
    if(exist('n','var'))  % get complete image
        x=zeros(n);
        x(maskIdx)=s;
    else
        if(min(size(maskIdx))>1)  % for maskFunc(s,mask);
            mask=maskIdx;
            if(length(s(:))==length(mask(:)))  % apply mask to image s;
                x=s(find(mask~=0));
            else  % restore the complete image
                x=zeros(size(mask));
                maskIdx=find(mask~=0);
                x(maskIdx)=s;
            end
        else  % apply the mask to image s;
            if(isempty(maskIdx))
                if(min(size(s))==1)
                    x=reshape(s,sqrt(length(s)),[]);
                else
                    x=s(:);
                end
            else
                x=s(maskIdx);
            end
        end
    end
end


