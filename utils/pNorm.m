function out=pNorm(x,p)
    switch (p)
        case 0
            out = (x~=0);
            out = sum(out(:));
        case 1
            out = sum(abs(x(:)));
        case 2
            out = norm(x,'fro');
        case inf
            out = max(abs(x(:)));
        otherwise
            out = abs(x).^p;
            out = sum(out(:))^(1/p);
    end
end
