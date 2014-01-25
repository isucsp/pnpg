%   Test the mex function parFwdPrj.c
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Wed 13 Nov 2013 06:36:14 PM CST

N=512;
img=zeros(N,N);
img(50:100,N/2-50:N/2+50)=1;
load(['../Mask' num2str(N) 'Circle.mat']);
addpath('../');
maskIdx=find(Mask);
img=img(maskIdx);
conf.bw=1; conf.nc=N; conf.nr=N; conf.prjWidth=N; conf.theta=0:179;
tic;
out=mParPrj(img,maskIdx-1,conf,'forward');
toc
figure(911); showImg(reshape(out,N,[]));
tic;
xxx=mParPrj(out(:),maskIdx-1,conf,'backward');
toc
figure(912); showImgMask(xxx,Mask);

% calculate the diagnal elements of the PhiPhi^t matrix
sinogram=zeros(N,1);

if(1)
diagnal=zeros(N,180);
diagnal1=zeros(N-1,180);
for j=0:179
    j
    conf.theta=j;
    for i=1:N
        sinogram=sinogram*0; sinogram(i)=1;
        out=mParPrj(sinogram,maskIdx-1,conf,'backward');
        %figure(913); showImgMask(out,Mask);
        diagnal(i,j+1)=out(:)'*out(:);
        if(i==1) outOld=out; continue; end
        diagnal1(i-1,j+1)=outOld(:)'*out(:);
        outOld=out;
    end
end
end
        
% calculate the correlation between X-rays from different angles
sinogram=sinogram*0; sinogram(N/2)=1;
intersection=zeros(179:179);
for i=0:178
    i
    conf.theta=i;
    out=mParPrj(sinogram,maskIdx-1,conf,'backward');
    for j=i+1:179
        conf.theta=j;
        out1=mParPrj(sinogram,maskIdx-1,conf,'backward');
        intersection(i+1,j)=out(:)'*out1(:);
        %figure(913); showImgMask(out1,Mask);
    end
end


