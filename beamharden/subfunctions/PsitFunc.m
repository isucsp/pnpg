function s=PsitFunc(x,Wt,m,maskIdx,wvltIdx)

% s: sparse vector which is taken from Wavelet transform of image x
% m, n: size of image

% converting wavelet coefficients (sparse representation) into samples
if(nargin==5)
    temp=zeros(m^2,1);
    temp(maskIdx)=x(:);
    x=reshape(temp,m,m);
else
    m=length(x(:));
    x=reshape(x,m,m);
end
s=Wt(x);
if(nargin==5)
    s=s(wvltIdx);
else
    s=s(:);
end
