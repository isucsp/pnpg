function maskk=wvltMask(Mask,dwt_L,daubcqfx,filename)
%imageName='realct';
%switch lower(imageName)
%    case 'phantom'
%        Img2D=phantom(256);
%        Mask=(Img2D~=0);
%        dwt_L=6;        %levels of wavelet transform
%        daubcqfx=2;
%    case 'brain'
%        Img2D = double(imread('brain.bmp'));
%        Mask=imread('brainMask.bmp');
%        dwt_L=8;        %levels of wavelet transform
%        daubcqfx=8;
%    case 'industry'
%        Img2D = double(imread('Industry1.bmp'));
%        Mask=imread('IndustryPbMask1L.bmp');
%        dwt_L=4;        %levels of wavelet transform
%        daubcqfx=8;
%    case 'realct'
%        Img2D = double(imread('RealCT.bmp'));
%        load('RealCTMask_03.mat');
%%       load 'Mask1024Circle.mat';
%        dwt_L=8;        %levels of wavelet transform
%        daubcqfx=6;
%    case 'simphantom'
%        load 'Phantom512.mat';
%        load 'Mask512Circle.mat';
%        dwt_L=5;
%        daubcqfx=2;
%end

Mask=(Mask~=0);
wav=daubcqf(daubcqfx);
W=@(z) midwt(z,wav,dwt_L);
Wt=@(z) mdwt(z,wav,dwt_L);
% coeff=Wt(Img2D);
% figure; semilogy(sort(abs(coeff(:)),'descend')); grid on;return;

test =zeros(size(Mask));
maskk=zeros(size(Mask));

maskIdx=find(Mask);
len=length(maskIdx);
str='';
tic; oldt=0;
for i=1:len
    newt=toc;
    if(abs(newt-oldt)>0.5) 
        strlen = length(str);
        str=sprintf('processing...%d/%d     ',i,len);
        fprintf([repmat('\b',1,strlen) '%s'],str);
        oldt=newt;
    end
    test(maskIdx(i))=1;
    maskk=maskk|Wt(test);
    test(maskIdx(i))=0;
end
fprintf('\nFinished\n');
if(~exist(filename,'file') && nargin>=4)
    save(filename,'maskk');
end
%save('RealCTMaskWvlt_03.mat','maskk');
%imwrite(maskk,'RealCTMaskWvlt_03.bmp','bmp');
%figure; showImg(maskk);
