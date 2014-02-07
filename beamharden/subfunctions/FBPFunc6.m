function x = FBPFunc6(y,theta,Ts,maskIdx)
% x: wavelet coefficients vector
% y: the input projection data in frequency domain
% f: coefficients of the filter applying on y
% st: 
% h: wavelet scaling filter
% L: level of wavelet decomposition
% m: size of image
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Fri 07 Feb 2014 02:46:21 AM CST

n=length(theta);
m=length(y(:))/n;
y=reshape(y,m,n);
y=[zeros(m/2,n); y; zeros(m/2,n)];
[Num_pixel,Num_proj]=size(y);

m_2D=[m, m];
J=[1,1]*3;                       % NUFFT interpolation neighborhood
K=2.^ceil(log2(m_2D*2));         % oversampling rate

r=pi*linspace(-1,1-2/Num_pixel,Num_pixel)';
xc=r*cos(theta(:)'*pi/180);
yc=r*sin(theta(:)'*pi/180);
om=[yc(:),xc(:)];
st=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');

f=designFilter('renliang1',m*2,Ts);

y=fftshift(fft(fftshift(y,1)),1)*Ts;
y=y.*repmat(f(:),1,Num_proj);
x=nufft_adj(y(:),st)*pi/Ts/st.M;

if(nargin==3)
    x=real(x(:));
else
    x=real(x(maskIdx));
end

