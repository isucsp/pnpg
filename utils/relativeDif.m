
function [dif,nm] = relativeDif(a,b,p)
    if(nargin==2) p=2; end
    nm=pNorm(a-b,p); normb=pNorm(b,p);
    dif = nm/max(normb,1e-10);
end

