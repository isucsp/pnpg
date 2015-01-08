function run_demo()

%Demonstration code for Poisson image deblurring using a Hessian Schatten-norm
%regularization approach. 

path=which('run_demo.m');
idx=strfind(path,filesep);
path=path(1:idx(end));
load([path 'demo_data']);

stream = RandStream('mcg16807', 'Seed',1230);
RandStream.setGlobalStream(stream);

%Create blurred and noisy version of the image
peak=20;
b=3; %bacground intensity.
fs=f/max(f(:))*peak; %scaling of the underlying image so as to simulate the 
% level of the poisson noise. The smaller the values of the original
% intensities the higher the noise level. 

y=imfilter(fs,h,'conv','circular')+b;
y=poissrnd(y);%Poisson noise. 


%Set of parameters for the deblurring algorithm (see HSPIRAL1.m and
% HSPIRAL2.m for a description.)

lambda=0.055; % Regularization parameter
options={'x_init',y,'alpha',10*lambda,'iter',200,'verbose',true,'showfig',...
  false,'tol',1e-5,'img',fs,'prox_iter',5,'bc','reflexive','bounds',[0 inf],'b',b};

%Deblurring using the HSPIRAL1 algorithm with the nuclear norm (Schatten 
%norm of order 1) Hessian Regularizer 
[x,fun_val,residual,ISNR]=HSPIRAL1(y,h,lambda,options{:},'snorm','nuclear');

figure(100);
imshow(fs,[]);title('Ground-truth','fontsize',16);
figure(101);
imshow(y,[]);title('Blurred and noisy image','fontsize',16);
figure(102);imshow(x,[]);title('Restored image','fontsize',16);


figure(103);plot(fun_val);ylabel('Objective function','fontsize',16);
xlabel('Number of iterations','fontsize',16);
title('Evolution of the objective cost funtion','fontsize',16);

figure(104);plot(residual);ylabel('ADMM Residual','fontsize',16);
xlabel('Number of iterations','fontsize',16);
title('Evolution of the residual in ADMM','fontsize',16);

figure(105);plot(ISNR);ylabel('ISNR','fontsize',16);xlabel('Number of iterations','fontsize',16);
title('Evolution of the SNR improvement','fontsize',16);set(gca,'fontsize',16)
