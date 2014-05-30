
R=300;
img = zeros(n,n);
xx=-n/2:n/2-1;
yy=xx;
[xx,yy]=meshgrid(xx,yy);

img(xx.^2+yy.^2<=R*R)=1;

conf.n=n; conf.prjWidth=n; conf.np=360; 
conf.prjFull=360; conf.dSize=1; conf.effectiveRate=1; 
conf.d=dist;

mPrj(0,conf,'config');
mPrj(0,0,'showConf');

Phi= @(s) mPrj(reshape(single(s),conf.n,[])',conf,'forward')*Ts;
sino=double(reshape(Phi(img),[],360))/Ts;

figure; plot(0:359,sino(n/2+1,:),'.');
