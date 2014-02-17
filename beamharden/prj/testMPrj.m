%   Test the mex function mPrj.c
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Sun 16 Feb 2014 05:10:32 PM CST

N=512;
img=zeros(N,N);
img(50:100,N/2-50:N/2+50)=1;

conf.n=N; conf.prjWidth=N;
conf.np=180; conf.prjFull=360;
conf.dSize=1.0; %(n-1)/(Num_pixel+1);
conf.effectiveRate=1.0; conf.d=0.0;

mPrj(0,conf,'config');

Phi =@(s) mPrj(s,0,'forward');
Phit=@(s) mPrj(s,0,'backward');

tic;
out=Phi(img);
toc
out=reshape(out,N,[]);
figure(911); showImg(out);
tic;
xxx=Phit(out);
toc
figure(912); showImg(xxx);

% test the symmetry of the matrices Phi and Phit
RPar=testTranspose(Phi,Phit,N*conf.np,N^2,'Phi',100);

return;

N=512;
img=zeros(N,N);
img(50:100,N/2-50:N/2+50)=1;

conf.n=N; conf.prjWidth=N;
conf.np=360; conf.prjFull=360;
conf.dSize=1.0; %(n-1)/(Num_pixel+1);
conf.effectiveRate=1.0; conf.d=10000.0;

mPrj(0,conf,'config');

Phi =@(s) mPrj(s,0,'forward');
Phit=@(s) mPrj(s,0,'backward');

tic;
out=Phi(img);
toc
out=reshape(out,N,[]);
figure(911); showImg(out);
tic;
xxx=Phit(out);
toc
figure(912); showImg(xxx);

% test the symmetry of the matrices Phi and Phit
RFan=testTranspose(Phi,Phit,N*conf.np,N^2,'Phi',100);

return;
% calculate the diagnal elements of the PhiPhi^t matrix
sinogram=zeros(N,180);

if(0)
diagnal=zeros(N,180);
diagnal1=zeros(N-1,180);
for j=0:179
    j
    conf.theta=j;
    for i=1:N
        sinogram=sinogram*0; sinogram(i)=1;
        out=Phit(sinogram,maskIdx-1,conf,'backward');
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


