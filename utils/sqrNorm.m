function out=sqrNorm(x)
    out = conj(x).*x;
    out = sum(out(:));
end
