function out = realInnerProd(a,b)
    out = real(conj(a).*b);
    out = sum(out(:));
end
