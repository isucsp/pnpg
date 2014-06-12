
function dif = relativeDif(a,b)
    normb=pNorm(b);
    if(normb==0) dif=0;
    else
        dif = pNorm(a-b)/normb;
    end
end

