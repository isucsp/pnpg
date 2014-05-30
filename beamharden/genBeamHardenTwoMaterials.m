%   Simulate and generate measurements for beamhardening effect
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Tue 16 Apr 2013 05:36:30 AM CDT

clear;
spark=1;
loadNIST
iron(:,1)=iron(:,1)*1e3;
A150(:,1)=A150(:,1)*1e3;
% from the original to 130 points, increase 19%
% from 130 points to 1300 points, increasing 1e-6, which is not necessary.
epsilon=20:1:150;
muIron=interp1(iron(10:end,1),iron(10:end,2),...
    epsilon,'spline');
muA150=interp1(A150(10:end,1),A150(10:end,2),...
    epsilon,'spline');
Ie=gampdf((epsilon-20)*16/100,5,1);
if(spark)
    Ie(45)=Ie(45)*1000;
    Ie(22)=Ie(22)*1000;
end

[temp,idx]=sort(muIron);
%figure; loglog(iron(:,1),iron(:,2));
%figure; semilogy(epsilon,muIron,'-');
%hold on;loglog(iron(:,1),iron(:,2),'*');
figure; loglog(iron(13:20,2),...
    gampdf((iron(13:20,1)-20)*16/100,5,1),'*');
hold on; semilogy(muIron,Ie);
%figure; semilogx(iron(:,1),Ie);
figure; semilogy(epsilon,Ie,'*-');

Ts=0.008;
Img2D=double(imread('binaryCasting.bmp'));
load('Mask1024Circle.mat');
maskIdx=find(Mask~=0);
conf.bw=1; conf.nc=1024; conf.nr=1024; conf.prjWidth=1024;
conf.theta=0:179;
temp=1:numel(Img2D);
Phi= @(s) mParPrj(s,temp-1,conf,'forward')*Ts;
Phit =@(s) mParPrj(s,temp-1,conf,'backward')*Ts;
PhiM=@(s) mParPrj(s,maskIdx-1,conf,'forward')*Ts;
PhiMt=@(s) mParPrj(s,maskIdx-1,conf,'backward')*Ts;
Num_proj=180;
theta_full=linspace(0,pi-pi/Num_proj,Num_proj); %The angle of projections
FBP=@(s) FBPFunc7(s,theta_full,1:Num_proj,Ts,maskIdx)*Ts;

img1=Img2D;
N=1024;
[c,r]=meshgrid(1:N);
c1=300;
r1=770;
rad=70;
img1((c-c1).^2+(r-r1).^2<rad^2)=2;

c1=600;
r1=845;
rad=70;
img1((c-c1).^2+(r-r1).^2<rad^2)=3;

img=img1; Img2D=img1;
save('twoMaterialCasting.mat','Img2D');
imgA150=zeros(size(img));
imgIron=imgA150;
imgA150(img==1)=densityA150;
imgIron(img>=2)=densityIron;
img(img==1)=densityA150;
img(img>=2)=densityIron;

epsilon=epsilon(:);
deltaEpsilon=mean([epsilon(:) [epsilon(2:end); epsilon(end)]],2)-...
    mean([epsilon(:) [epsilon(1); epsilon(1:end-1)]],2);
PhiAlphaA150=Phi(imgA150);
PhiAlphaIron=Phi(imgIron);

Imea=0;
for i=1:length(epsilon)
    Imea=Imea+exp(-PhiAlphaIron*muIron(i)-PhiAlphaA150*muA150(i))...
        *Ie(i)*deltaEpsilon(i);
end

y=-log(Imea/max(Imea(:)));
rec=FBP(y);
rec2=FBP(Phi(img));
CTdata_mtrx=reshape(Imea,1024,[]);
if(spark)
    save('CTdata_220TwoMaterialsSpark.mat','CTdata_mtrx');
else
    save('CTdata_220TwoMaterials.mat','CTdata_mtrx');
end

