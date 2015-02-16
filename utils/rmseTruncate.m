function rmse=rmseTruncate(x,trueAlpha)
    if(nargin<2) trueAlpha=x.opt.trueAlpha; end;
    alpha=x.alpha; alpha(alpha<0)=0;
    rmse=sqrNorm(alpha-trueAlpha)/sqrNorm(trueAlpha);
end

