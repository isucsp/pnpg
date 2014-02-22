function x = FBPFunc8(y,conf,Ts,maskIdx)
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
%   $Revision: 0.1 $ $Date: Mon 17 Feb 2014 10:35:42 PM CST

n=conf.np;
m=length(y(:))/n;
y=reshape(y,m,n);
[Num_pixel,Num_proj]=size(y);

temp=2^(ceil(log2(Num_pixel)));
y=[zeros(ceil((temp*2-Num_pixel)/2),Num_proj); y;...
    zeros(floor((temp*2-Num_pixel)/2),Num_proj)];
[Num_pixel,Num_proj]=size(y);

f=designFilter('renliang1',Num_pixel,Ts);

y=fftshift(fft(fftshift(y,1)),1)*Ts;
y=y.*repmat(f(:),1,Num_proj);
y=fftshift(ifft(fftshift(y,1)),1)/Ts;
y=real(y);
y=y(Num_pixel/2-floor(m/2):Num_pixel/2+floor((m-1)/2),:);
[Num_pixel,Num_proj]=size(y);

x=reshape(mPrj(y,0,'backward'),conf.n,[]);

