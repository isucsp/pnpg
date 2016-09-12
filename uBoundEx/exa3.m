clear all; close all;
% L(x)=0.5*||x||_2^2
% C={x| ||x-[2;0]||_2^2<=2 }

xstar=ones(2,1);
Psi=@(x)[x;-x];
Psit=@(x) x(1)-x(2);
%Ncx=[-a,a] with a>=0
Pncx=@(x) min(0,(x(1)-x(2))/2)*[1;-1];
g=[1;1];

uTrue=Inf

u=uBound(Psi,Psit,Pncx,xstar,g);

