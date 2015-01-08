function out=pNorm(x,p)
    if(nargin==1) p=2; end
    switch (p)
        case 0
            out = (p~=0);
            out = sum(out(:));
        case 1
            out = sum(abs(x(:)));
        case 2
            out = x.*x;
            out = sqrt(sum(out(:)));
        case inf
            out = max(abs(x(:)));
        otherwise
            out = abs(x).^p;
            out = sum(out(:))^(1/p);
    end
end
