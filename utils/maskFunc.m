
function x = maskFunc(s,maskIdx,n)
    if(nargin==3)
        x=zeros(n,n);
        x(maskIdx)=s;
    else
        x=s(maskIdx);
    end
end


