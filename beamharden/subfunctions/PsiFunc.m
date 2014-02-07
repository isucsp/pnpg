function x=PsiFunc(s,W,m,maskIdx,wvltIdx)

% s: sparse vector which is taken from Wavelet transform of image x
% m, n: size of image

% converting wavelet coefficients (sparse representation) into samples
if(nargin==5)
    temp=zeros(m^2,1);
    temp(wvltIdx)=s(:);
    s=reshape(temp,m,m);
else
    m=length(s(:));
    s=reshape(s,m,m);
end
x=W(s);
if(nargin==5)
    x=x(maskIdx);
else
    x=x(:);
end

