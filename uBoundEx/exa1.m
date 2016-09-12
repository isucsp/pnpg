clear all; close all;
% L(x)=x_1+I_+(x_1)
% Psi = @(x) x;
% C={x| ||x-1||<=1 }

xstar=(1-sqrt(2)/2)*ones(2,1);
Psi=@(x)x;
Psit=Psi;
%Ncx=[-a,-a] with a>=0
Pncx=@(x) min(0,mean(x))*ones(2,1);
g=[1;0];

uTrue=Inf

u=uBound(Psi,Psit,Pncx,xstar,g);

