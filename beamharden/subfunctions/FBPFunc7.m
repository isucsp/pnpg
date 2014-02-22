function x = FBPFunc7(y,prjFull,prjNum,Ts,maskIdx)
% FBPFunc7      Filtered back projection reconstruction
%   This routine is implemented by the external C language. The C code 
%   implement the projection and back projection by the discrete methods.
%   x: wavelet coefficients vector
%   y: the input projection data in frequency domain
%   f: coefficients of the filter applying on y
%   st: 
%   h: wavelet scaling filter
%   L: level of wavelet decomposition
%   m: size of image
%   out = beamharden(***)
%   Phi         The projection matrix implementation function handle
%   Phit        Transpose of Phi
%   Psi         Inverse wavelet transform matrix from wavelet coefficients
%               to image.
%   Psit        Transpose of Psi
%   y           Log scale of Beamhardening measurement y=-log(I^{mea}/I_0)
%   xInit       Initial value for the algorithm
%   opt         Structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Mon 17 Feb 2014 08:05:42 PM CST

n=length(prjNum);
m=length(y(:))/n;
y=reshape(y,m,n);
y=[zeros(m/2,n); y; zeros(m/2,n)];
[Num_pixel,Num_proj]=size(y);

f=designFilter('renliang1',m*2,Ts);

y=fftshift(fft(fftshift(y,1)),1)*Ts;
y=y.*repmat(f(:),1,Num_proj);

y=fftshift(ifft(fftshift(y,1)),1)/Ts;
y=real(y);
y=y(m/2+1:m/2+m,:);
[Num_pixel,Num_proj]=size(y);
conf.bw=1; conf.nc=Num_pixel; conf.nr=Num_pixel; conf.prjWidth=Num_pixel;
conf.theta=(0:prjNum-1)*360/prjFull;
out=mParPrj(y,maskIdx-1,conf,'backward');
x=zeros(m,m);
x(maskIdx)=out(1:length(maskIdx));

