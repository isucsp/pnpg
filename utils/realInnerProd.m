function out = realInnerProd(a,b)
    out = real(a).*real(b)+imag(a).*imag(b);
    out = sum(out(:));
end
