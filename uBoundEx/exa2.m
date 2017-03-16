clear all; close all;
% L(x)=||x||_2^2
% Psi = @(x) x;
% C={x| ||x-1||<=1 }

xstar=(1-sqrt(2)/2)*ones(2,1);
Psi=@(x)x;
Psit=Psi;
%Ncx=[-a,-a] with a>=0
Pncx=@(x) min(0,mean(x))*ones(2,1);
g=(2-sqrt(2))*ones(2,1);

u=uBound(Psi,Psit,[],Pncx,xstar,g)

Psi=@(x)[x;-x];
Psit=@(x) x(1)-x(2);

u=uBound(Psi,Psit,[],Pncx,xstar,g)


