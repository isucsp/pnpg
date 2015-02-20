function rmse=rmseTruncate(x,trueAlpha)
    if(nargin<2) trueAlpha=x.opt.trueAlpha; end;
    if(isstruct(x))
        alpha=x.alpha;
    else
        alpha=x(:);
    end
    alpha(alpha<0)=0;
    rmse=sqrNorm(alpha-trueAlpha)/sqrNorm(trueAlpha);
end

