clear;

Ts=0.1;
down=1;
sg = sino_geom('par', 'nb', 511, 'na', 720, 'dr',Ts, ...
    'orbit', 180, ... % curiously, 360 has less aliasing artifact
    'offset_r', 0, 'orbit_start', 0, 'down', down);
ig = image_geom('nx', 512, 'ny', 512, 'dx',Ts,'down', down);
ig.dy = ig.dx; % DF needs this for now

f.erot = -90;
ell=[...
    0       0       0.92	0.69	90	1;
    0       -0.0184	0.874	0.6624	90	-0.8;
    0.22	0       0.31	0.11	72	-0.2;
    -0.22	0       0.41	0.16	108	-0.2;
    0       0.35	0.25	0.21	90	0.1;
    0       0.1     0.046	0.046	0	0.1;
    0       -0.1	0.046	0.046	0	0.1;
    -0.08	-0.605	0.046	0.023	0	0.1;
    0       -0.605	0.023	0.023	0	0.1;
    0.06	-0.605	0.046	0.023	90	0.1];

ell(:,1:4) = ell(:,1:4) * ig.fov/2;
% ell = [0 0 200 200 0 0; 13 13 10 10 0 1000];
%figure; imshow(xtrue,[]);
CTdata = ellipse_sino(sg, ell, 'oversample',1);
CTdata=[CTdata(:,1), CTdata(end:-1:1,end:-1:2)]/Ts;
ell(:,1)=ell(:,1)+Ts/2; ell(:,2)=ell(:,2)-Ts/2;
[xtrue ell1] = ellipse_im(ig, ell, 'oversample', 1, 'rot', f.erot);

Nratio=[]; i=0;
% choose which mask to use
imageName='Phantom';  % also can be 'phantom' and 'brain'

%Select projections
spacing=1; %The spacing of projections
% Restricted angles
%theta_idx=[5:40,45:120,125:160,165:180];
theta_idx=1:2:720;
i=i+1;

%filter design
filter='renliang'; %'shepp-logan'; %'ram-lak'; %

saveImg=1;

switch lower(imageName)
    case 'phantom'
        Img2D=xtrue;
        Mask=(Img2D~=0);
        load('wvltMaskPhantom.mat');
        dwt_L=6;        %levels of wavelet transform
        wav=daubcqf(2);
    case 'industry'
        load 'IndustryPb1.mat'
        wav=daubcqf(8);
        dwt_L=4;
        Img2D = double(imread('Industry1.bmp'));
        fprintf('Loading Mask ...\n');
        Mask=imread('IndustryPbMask1L.bmp');
        load('IndustryPbMaskWvlt1L.mat');
end

[N,M]=size(CTdata);
CTdata=CTdata((N+1)/2-255:(N+1)/2+255,:);

theta_full=linspace(0,pi-pi/M,M); %The angle of projections
theta=theta_full(theta_idx);
CTdata=CTdata(:,theta_idx);
[Num_pixel,Num_proj]=size(CTdata);

N_freq=2^ceil(log2(Num_pixel));
CTdata=[zeros(N_freq/2-(Num_pixel-1)/2,Num_proj); CTdata;...
    zeros(N_freq/2-(Num_pixel+1)/2,Num_proj)]/sqrt(N_freq);
y1=CTdata(:);
N=length(y1);
Num_pixel=N_freq;
m=Num_pixel^2;
mx=Num_pixel; my=Num_pixel;
% N_freq=768;
% N_freq=ceil(sqrt(mx^2+my^2)/2)*2;
%N_freq=Num_pixel*2;

% Zero freq at f_coeff(0)
% This is a point to exploit for better performance.
f_coeff=sqrt(designFilter(filter,Num_pixel)*pi/Num_proj*N_freq);

y2=[zeros(Num_pixel/2,Num_proj); CTdata; ...
    zeros(Num_pixel/2,Num_proj)];
y2=fft(fftshift(y2,1))/sqrt(2*Num_pixel);
y2=y2.*repmat(f_coeff,1,Num_proj);
y2=fftshift(ifft(y2),1)*sqrt(2*Num_pixel);
y2=y2(N_freq/4+1:3*N_freq/4,:);
y2=real(y2(:));

m_2D=[mx, my];
J=[1,1]*4;                       % NUFFT interpolation neighborhood
K=2.^ceil(log2(m_2D*2));         % oversampling rate

r=pi*linspace(-1,1-2/N_freq,N_freq)';
%The locations of measurements in Frequency plane
xc=r*cos(theta(:)');
yc=r*sin(theta(:)');
om=[yc(:),xc(:)];
stAdj=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');

r=pi*linspace(-1,1-2/Num_pixel,Num_pixel)';
%The locations of measurements in Frequency plane
xc=r*cos(theta(:)');
yc=r*sin(theta(:)');
om=[yc(:),xc(:)];
stFwd=nufft_init(om,m_2D,J,K,m_2D/2,'minmax:kb');

Phi1=@(s) PhiFunc1(s,f_coeff,stFwd,Num_pixel);
Phi2=@(s) PhiFunc2(s,f_coeff,stFwd,Num_pixel);

Phit=@(s) FBPFunc(s,f_coeff.^2,stAdj,Num_pixel,N_freq);
Phit1=@(s) PhitFunc1(s,f_coeff,stAdj,Num_pixel,N_freq);
Phit2=@(s) PhitFunc2(s,f_coeff,stAdj,Num_pixel,N_freq);

return;

yy1=Phi1(Img2D);
yy2=Phi2(Img2D);
xx1=Phit1(y1);
xx1=reshape(xx1,my,mx);
xx2=Phit2(y2);
xx2=reshape(xx2,my,mx);
figure; subplot(2,1,1);
plot(y1(1:512)); hold on; plot(yy1(1:512),'r');
subplot(2,1,2);
plot(y2(1:512)); hold on; plot(yy2(1:512),'r');

figure; subplot(2,1,1);
plot(Img2D(200,:)); hold on; plot(xx1(200,:),'r');
subplot(2,1,2);
plot(Img2D(200,:)); hold on; plot(xx2(200,:),'r');

return;

j=0;
prefix='BackProj';
j=j+1;
fprintf('%s, i=%d, j=%d\n',prefix,i,j);
s_BackProj=Phit(y);
Phis=Phi(s_BackProj);
res_BackProj(i,j)=norm(y-Phis);
Phis=reshape(Phis,[Num_pixel,Num_proj]);
Img2D_BackProj=reshape(s_BackProj,[my mx]);
error=norm(Img2D(:)-Img2D_BackProj(:));
figure;
subplot(2,2,1); showImg(Img2D_BackProj);
subplot(2,2,3); showImg((Img2D_BackProj-Img2D));
subplot(2,2,2); showImg(Phis);
subplot(2,2,4); showImg((Phis-CTdata));
figure; subplot(2,1,1);
plot(Img2D_BackProj(257,:),'r'); hold on;
plot(Img2D(257,:)); ylim([-0.2 1.2]);

subplot(2,1,2); plot(Img2D_BackProj(208,:),'r'); hold on;
plot(Img2D(208,:)); ylim([-0.2 1.2]);
