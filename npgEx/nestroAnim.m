
close all;
clear all;

A=[1,0.21;0.5,3];
x0=[0.1;0.1];
y=poisson(A*x0*100)/100;
y=A*x0;
f=@(a,b) (A(1,1)*a+A(1,2)*b-y(1)+A(2,1)*a+A(2,2)*b-y(2))...
    +y(1)*log(y(1)./(A(1,1)*a+A(1,2)*b))...
    +y(2)*log(y(2)./(A(2,1)*a+A(2,2)*b));
[x1,x2]=meshgrid(linspace(0,3*x0(1),1000), linspace(0,3*x0(2),1000));
fx=f(x1,x2);
Phi=@(x) A*x;
Phit=@(x) A.'*x;
Psi=@(x) x;

%xInit=Phit(y);
xInit=[0.1,1];
xInit=[6,5];
xInit=[0,3];
xInit=[1,0];
xInit=[5.5,0.5];
xInit=[4,5];
xInit=[3,0];
xInit=[1,1]*1e-5;

opt.maxItr=5e4; opt.thresh=1e-6; opt.debugLevel=1;
opt.u=1e-19;
opt.L = max(svd(A))^2;
opt.continuation=false; opt.alphaStep='pnpg'; opt.saveXtrace=true;
opt.noiseType='Poisson';
opt.adaptiveStep=true;
out=solver(Phi,Phit,Psi,Psi,y,xInit,opt);
out.alpha

v=logspace(-2,4,20);
figure; contour(x1,x2,fx,v);
hold on; plot(xInit(1),xInit(2),'r.');

for i=1:out.p
    plot(out.alphaTrace(1,i),out.alphaTrace(2,i),'r.');
    pause(0.2);
end


keyboard

opt.alphaStep='pnpg'; opt.adaptiveStep=false;
out=solver(Phi,Phit,Psi,Psi,y,xInit,opt);
out.alpha
figure; contour(x1,x2,fx,v);
hold on; plot(xInit(1),xInit(2),'r.');
for i=1:out.p
    plot(out.alphaTrace(1,i),out.alphaTrace(2,i),'r.');
    %pause(0.2);
end

keyboard

opt.alphaStep='pg'; opt.adaptiveStep=true;
out=solver(Phi,Phit,Psi,Psi,y,xInit,opt);
out.alpha
figure; contour(x1,x2,fx,v);
hold on; plot(xInit(1),xInit(2),'r.');
for i=1:out.p
    plot(out.alphaTrace(1,i),out.alphaTrace(2,i),'r.');
    %pause(0.2);
end



