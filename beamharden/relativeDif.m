
function [dif,nm] = relativeDif(a,b,p)
    if(nargin==2) p=2; end
    nm=pNorm(a-b,p); normb=pNorm(b,p);
    if(nm==0 && normb==0)
        dif=0;
    else
        dif = nm/normb;
    end
end

