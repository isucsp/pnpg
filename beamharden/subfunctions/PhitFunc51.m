function x = PhitFunc51(y,f,st,m,Ts,maskIdx)
% x: wavelet coefficients vector
% y: the input projection data in frequency domain
% f: coefficients of the filter applying on y
% st: 
% h: wavelet scaling filter
% L: level of wavelet decomposition
% m: size of image

% Ts works to keep the right scale and unit

y=reshape(y,st.Num_pixel,st.Num_proj);
y=fftshift(fft(fftshift(y,1)),1)/st.Num_pixel/Ts;
x=nufft_adj(y(:),st)*Ts^2;
if(nargin>5)
    x=real(x(maskIdx));
end
x=real(x(:));

