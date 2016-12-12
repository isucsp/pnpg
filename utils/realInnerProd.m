function out = realInnerProd(a,b)
    if(isreal(a) && isreal(b))
        out = a.*b;
    else
        out = real(a).*real(b)+imag(a).*imag(b);
    end
    out = sum(out(:));
end
