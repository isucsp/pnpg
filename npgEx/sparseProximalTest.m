clear all; close all;

[~,~,~,Psi,Psit,OPT,~,~]=loadLinear([],0);

maxItr=4e3; thresh=1e-10;
u=1;
pOut=[];
prj_C=@(x) max(x,0);

for i=1:5

    a=randn(size(OPT.trueX))*10^(i-3)+OPT.trueX;

    opt=[];
    opt.outLevel=1;
    opt.debugLevel=2;
    opt.Lip=@(u)u^2;
    opt.initStep='fixed';
    opt.adaptiveStep=false;
    opt.backtracking=false;
    opt.dualGap=true;
    tStart=tic;
    sp1=sparseProximal(Psi,Psit,prj_C,'pnpg',opt);
    for j=1:1
        init=[];
        [x1,itr1,pOut,out1]=sp1.prox(a,u,thresh,maxItr,init);
    end
    t1=toc(tStart);

    tStart=tic;
    for j=1:1
        init=[];
        [x2,itr2,pOut] = admm(Psi,Psit,a,u,thresh,maxItr,false,init);
    end
    t2=toc(tStart);

    opt.debugLevel=0;
    opt.debugLevel=1;
    opt.debugLevel=2;
    opt=[];
    tStart=tic;
    sp3=sparseProximal(Psi,Psit,prj_C,'admm',opt);
    for j=1:1
        init=[];
        [x3,itr3,pOut3]=sp3.prox(a,u,thresh,maxItr,init);
    end
    t3=toc(tStart);


    f=@(x)0.5*sqrNorm(x-a)+u*sp3.val(x);
    fprintf('objective=%14.10g, itr=%d, gap=%10.8g\n',f(x1),itr1,out1.gap);
    fprintf('objective=%14.10g, itr=%d\n',f(x2),itr2);
    fprintf('objective=%14.10g, itr=%d\n',f(x3),itr3);

    fprintf('results: t1=%g, t2=%g, t3=%g\n', t1, t2, t3);

%    keyboard
end

