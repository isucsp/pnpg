function x = FBPFunc6(y,theta_full,theta_idx,Ts,maskIdx)
% x: wavelet coefficients vector
% y: the input projection data in frequency domain
% f: coefficients of the filter applying on y
% st: 
% h: wavelet scaling filter
% L: level of wavelet decomposition
% m: size of image
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Sun Sep 16 19:15:24 CDT 2012

n=length(theta_idx(:));
m=length(y(:))/n;
y=reshape(y,m,n);
y=[zeros(m/2,n); y; zeros(m/2,n)];
[Num_pixel,Num_proj]=size(y);

m_2D=[m, m];
J=[1,1]*3;                       % NUFFT interpolation neighborhood
K=2.^ceil(log2(m_2D*2));         % oversampling rate

theta=theta_full(theta_idx);
r=pi*linspace(-1,1-2/Num_pixel,Num_pixel)';
xc=r*cos(theta(:)');
yc=r*sin(theta(:)');
om=[yc(:),xc(:)];
st=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');

f=designFilter('renliang1',m*2,Ts);

y=fftshift(fft(fftshift(y,1)),1)*Ts;
y=y.*repmat(f(:),1,Num_proj);
x=nufft_adj(y(:),st)*pi/Ts/st.M;

if(nargin==4)
    x=real(x(:));
else
    x=real(x(maskIdx));
end
