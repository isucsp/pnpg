function [x,fun_val,residual,ISNR]=HSPIRAL2(y,h,lambda,varargin)

%Hessian Schatten norm Poisson Image Reconstruction with Augmented
%Lagrangian (HSPIRAL)
%
% Forward model: y=P(Ax+b).  
% P: Poisson noise, A: circulant convolution, b: background constant.
%
% ========================== INPUT PARAMETERS (required) ==================
% Parameters     Values description
% =========================================================================
% y             Noisy blured image.
% h             Blur kernel (PSF).
% lambda        Regularization penalty.
% ======================== OPTIONAL INPUT PARAMETERS ======================
% Parameters    Values description
%
% img           Original Image. (For the compution of the ISNR improvement)
% alpha         Quadratic penalty parameter (Default: lambda)
% x_init        Initial guess for x. (Default: [])
% iter          Number of iterations (Default: 100)
% tol           Stopping threshold: Relative normed difference between two
%               successive iterations (Default:1e-5)
% verbose       If verbose is set on then info for each iteration is
%               printed on screen. (Default: false)
% showfig       If showfig is set on the result of the deconvolution in
%               each iteration is shown on screen. (Default: false)
% snorm         Specifies the type of the Hessian Schatten norm.
%               {'spectral'|'nuclear'|'frobenius'|'Sp'}. (Default:
%               'frobenious').
% order         In case of a general Sp-norm regularizer it specifies the
%               order of the norm. (Default: [])
% bounds        Minimize the Objective involving the Hessian Schatten Norm
%               with a box constraint on the solution (Default: [0 +inf])
% b             Constant value for the background (Default: 0)
% =========================================================================
% ========================== OUTPUT PARAMETERS ============================
% x             Deconvolved image.
% fun_val       The evolution of the objective function.
% residual      The evolution of the residual in the ADMM iterations. 
% ISNR          The evolution of the SNR increase (if the original image is
%               provided.)
% =========================================================================
%
% Author: stamatis.lefkimmiatis@epfl.ch
%
% =========================================================================

%%Example
% peak=10;
% x=double(imread('cameraman.tif'));
% x=x./max(x(:))*peak;
% h=ones(9)/81;
% y=imfilter(x,h,'conv','circular');
% y=poissrnd(y.*mask);
%[xest,fun_val,residual,ISNR]=HSPIRAL2(y,h,0.025,'alpha',.1,'verbose',...
%true,'x_init',y,'iter',200,'verbose',true,'img',x,'snorm','frobenius','mask',mask);


[x_init,alpha,iter,verbose,showfig,tol,img,snorm,order,bounds,b]=...
  process_options(varargin,'x_init',[],'alpha',lambda,'iter',100,'verbose',...
  false,'showfig',false,'tol',1e-5,'img',[],'snorm','frobenius','order',...
  [],'bounds',[0 +inf],'b',0);

if ~(isequal(snorm,'frobenius') || isequal(snorm,'spectral') ...
    || isequal(snorm,'nuclear') || isequal(snorm,'Sp'))
  error('HSPIRAL2:: Unknown type of norm.');
end

if isequal(snorm,'Sp') && isempty(order)
  error('HSPIRAL2:: You have to specify the order of the Sp-norm.');
end

if ~isempty(order) && order < 1
  error('HSPIRAL2:: The order of the Sp-norm must be >=1.');
end

if b < 0
  error('HSPIRAL2: Background must have a non-negative value');
end


% Initializations
hdim=size(h);
ydim=size(y);


%zero padding h to match the image size.
hzp=padarray(h,ydim-hdim,'post');

%Circular shift of the PSF so as to have consistent results by applying
%either one of the following 2 operations on an image x (hzp is the zero
%padded psf)
%1)imfilter(x,h,'conv','circular')
%2)ifft2(fft2(x).*fft2(hzp));

%============================ IMPORTANT NOTE ==============================
%Matlab's imfilter picks as the origin of h its central value, i.e
% chy=floor((hr+1)/2) and chx=floor((hc+1)/2). Therefore without circular
% shifting hzp, the above 2 operations are not consistent.
%==========================================================================

hzp=circshift(hzp,-floor(hdim/2));
H=fft2(hzp);
clear hzp h

%x0 Initialization
Bop=@(y)real(ifft2(H.*fft2(y)));%Function handle that corresponds to
%the blurring operation of y.

%H*H=Dxx*Dxx+2Dxy*Dxy+Dyy*Dyy, which is a circulant operator since each one
%of the suboperators are circulant. The addition and multiplication of
%circulant operators results to a circulant operator.

Dxx=[1 -2 1 0 0]';% Dxx Operator
Dxx=padarray(Dxx,[0 2]);
Dyy=Dxx'; %Dyy Operator
Dxy=[1 -1 0;-1 1 0;0 0 0];%Dxy Operator
Dxy=padarray(Dxy,[1 1]);
center=[3 3];

Dxx(ydim(1),ydim(2))=0;Dxy(ydim(1),ydim(2))=0;Dyy(ydim(1),ydim(2))=0;

%Circular shift of the operator T1 so as to have consistent results by
%applying either one of the following 2 operations on an image x (T1zp is
%the zero padded operator)
%1)imfilter(x,T1,'conv','circular')
%2)ifft2(fft2(x).*fft2(T1zp));

Dxx=fft2(circshift(Dxx,1-center));
Dxy=fft2(circshift(Dxy,1-center));
Dyy=fft2(circshift(Dyy,1-center));

%K=H'*H+Hop'*Hop+I in the Fourier domain where Hop the Hessian Operator.
K=abs(H).^2+abs(Dxx).^2+2*abs(Dxy).^2+abs(Dyy).^2+1;

s1=zeros(ydim);%Lagrangian multipliers.
s2=zeros([ydim,3]);%Lagrangian multipliers.
s3=zeros(ydim);%Lagrangian multipliers.

if isempty(x_init)
  x_init=y;
end
x=x_init;

fun_val=zeros(iter,1);
ISNR=zeros(iter,1);
residual=zeros(iter,1);
bc='circular';

count=0;
if verbose
  fprintf('\t\t***********************************\n');
  fprintf('\t\t**  Deconvolution with ADMM      **\n');
  fprintf('\t\t***********************************\n');
  fprintf('#iter       fun-val      relative-dif     residual       ISNR\n')
  fprintf('================================================================\n');
end
for i=1:iter
  
  %u1=Bop(x)+b+s1;
  %u2=Hessianop2D(x)+s2;
  %u3=x+s3;
  z1=Poissonshrink(y,Bop(x)+b+s1,alpha);
  z2=Hessianshrink(HessianOp2D(x,bc)+s2,lambda/alpha,snorm,order);
  z3=project(x+s3,bounds);
  
  %v1=z1-b-s1;
  %v2=z2-s2;
  %v3=z3-s3;
  R=fft2(z2-s2);
  xnew=real(ifft2((conj(H).*fft2(z1-b-s1)+conj(Dxx).*R(:,:,1)+...
    2*conj(Dxy).*R(:,:,2)+conj(Dyy).*R(:,:,3)+fft2(z3-s3))./K));
  %xnew(xnew<0)=0;
  
  s1=s1+Bop(xnew)+b-z1;
  s2=s2+HessianOp2D(xnew,bc)-z2;
  s3=s3+xnew-z3;
  
  fun_val(i)=cost(y,Bop,xnew,lambda,snorm,order,bc,b);
  
  re=norm(xnew-x,'fro')/norm(x,'fro');%relative error
  if (re < tol)
    count=count+1;
  else
    count=0;
  end
  
  if verbose
    if ~isempty(img)
      ISNR(i)=20*log10(norm(y-img,'fro')/norm(xnew-img,'fro'));
      temp=HessianOp2D(xnew,bc)-z2;
      temp=temp(:,:,1).^2+2*temp(:,:,2).^2+temp(:,:,3).^2;
      temp=sum(temp(:));
      residual(i)=sqrt(norm(Bop(xnew)+b-z1,'fro')^2+temp+norm(xnew-z3,'fro')^2);
      % printing the information of the current iteration
      fprintf('%3d \t %3.5f \t %3.5f \t %3.5f \t %3.5f\n',i,fun_val(i),re,residual(i),ISNR(i));
    else
      temp=HessianOp2D(xnew,bc)-z2;
      temp=temp(:,:,1).^2+2*temp(:,:,2).^2+temp(:,:,3).^2;
      temp=sum(temp(:));
      residual(i)=sqrt(norm(Bop(xnew)+b-z1,'fro')^2+temp+norm(xnew-z3,'fro')^2);
      fprintf('%3d \t %3.5f \t %3.5f \t %3.5f\n',i,fun_val(i),re,residual(i));
      %fprintf('%3d \t %3.5f\t%1.5f\n',i,fun_val(i),re);
    end
  end
  x=xnew;
    
  if showfig
    fh=figure(1);
    figure(fh);
    msg=['iteration: ' num2str(i) ' ,ISNR: ' num2str(ISNR(i))];
    set(fh,'name',msg);imshow(x,[]);
  end
  
  if count >=5
    ISNR(i+1:end)=[];
    fun_val(i+1:end)=[];
    residual(i+1:end)=[];
    break;
  end
  
end

function [Q,Hnorm]=cost(y,Bop,f,lambda,snorm,p,bc,b)

fxx=(f-2*shift(f,[-1,0],bc)+shift(f,[-2,0],bc));
fyy=(f-2*shift(f,[0,-1],bc)+shift(f,[0,-2],bc));
fxy=(f-shift(f,[0,-1],bc)-shift(f,[-1,0],bc)+shift(f,[-1,-1],bc));

switch snorm
  case 'spectral'
    %Sum of the Hessian spectral radius
    Lf=fxx+fyy;%Laplacian of the image f
    Of=sqrt((fxx-fyy).^2+4*fxy.^2);% Amplitude of the Orientation vector
    Hnorm=sum(sum(0.5*(abs(Lf)+Of)));%Sum of the Hessian spectral radius
  case 'frobenius'
    Hnorm=sum(sqrt(fxx(:).^2+fyy(:).^2+2*fxy(:).^2));
  case 'nuclear'
    Lf=fxx+fyy;%Laplacian of the image f
    Of=sqrt((fxx-fyy).^2+4*fxy.^2);% Amplitude of the Orientation vector
    Hnorm=0.5*sum(abs(Lf(:)+Of(:))+abs(Lf(:)-Of(:)));
  case 'Sp'
    e1=0.5*(fxx+fyy+sqrt((fxx-fyy).^2+4*fxy.^2));
    e2=0.5*(fxx+fyy-sqrt((fxx-fyy).^2+4*fxy.^2));
    Hnorm=sum((abs(e1(:)).^p+abs(e2(:)).^p).^(1/p));
  otherwise
    error('HSPIRAL2::Unknown type of norm.');
end

Hf=Bop(f)+b;
Hf(Hf<0)=0;
data_term=Hf(:)-y(:).*log(Hf(:));%0*log(0)=0 by convention. Matlab returns
%NAN so we have to fix this.
data_term(isnan(data_term))=0;
data_term(isinf(data_term))=0;
Q=sum(data_term)+lambda*Hnorm;


function z=Poissonshrink(y,u,alpha)
%It solves component-wise the minimization problem
%
% argmin {0.5*(z_i-u_i)^2+ (1/alpha)*(z_i-y_i*log(z_i))
%   z>=0
%------------------------------------------------------

u=u-(1/alpha);
z=.5*(u+sqrt(u.^2+(4*y/alpha)));


function X=Hessianshrink(Y,alpha,snorm,order)

switch snorm
  case 'frobenius'
    X=proxSpMat2x2(Y,2,alpha);
  case 'nuclear'
    X=proxSpMat2x2(Y,1,alpha);
  case 'spectral'
    X=proxSpMat2x2(Y,inf,alpha);
  case 'Sp'
    X=proxSpMat2x2(Y,order,alpha,1);
end

function x=project(y,bounds) %Projection onto the non-negative set
x=y;
x(x<bounds(1))=bounds(1);
x(x>bounds(2))=bounds(2);


function Hf=HessianOp2D(f,bc)

[r,c]=size(f);

fxx=(f-2*shift(f,[-1,0],bc)+shift(f,[-2,0],bc));
fyy=(f-2*shift(f,[0,-1],bc)+shift(f,[0,-2],bc));
fxy=(f-shift(f,[0,-1],bc)-shift(f,[-1,0],bc)+shift(f,[-1,-1],bc));

%Compute Hf (Apply to image f the Hessian Operator)
%Hf will be a cube (3D image)
Hf=zeros(r,c,3);
Hf(:,:,1)=fxx;
Hf(:,:,2)=fxy;
Hf(:,:,3)=fyy;
