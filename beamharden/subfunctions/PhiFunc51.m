function y = PhiFunc51(x,f,st,m,Ts,mask)
% y: the output measurements
% x: the input wavelet coefficients vector
% f: coefficients of the filter applying on y
% st: 
% h: wavelet scaling filter
% L: level of wavelet decomposition
% m: size of image
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Thu 07 May 2015 10:27:04 AM CDT

if(nargin>5)
    u=mask.b(x);
else
    u=reshape(x,m,m);
end
y=nufft(u,st)*Ts^2;
y=reshape(y,st.Num_pixel,st.Num_proj);
y=fftshift(y,1);
y=fftshift(ifft(y),1)/Ts;
y=real(y(:));

