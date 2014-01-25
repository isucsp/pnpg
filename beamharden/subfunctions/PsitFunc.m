function s=PsitFunc(x,Wt,m,maskIdx,wvltIdx)

% s: sparse vector which is taken from Wavelet transform of image x
% m, n: size of image

% converting wavelet coefficients (sparse representation) into samples
if(nargin==5)
    temp=zeros(m,1);
    temp(maskIdx)=x(:);
    x=reshape(temp,sqrt(m),sqrt(m));
else
    m=length(x(:));
    x=reshape(x,sqrt(m),sqrt(m));
end
s=Wt(x);
if(nargin==5)
    s=s(wvltIdx);
else
    s=s(:);
end
