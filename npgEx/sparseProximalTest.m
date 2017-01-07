clear all; close all;

[~,~,~,Psi,Psit,OPT,~,~]=loadLinear([],0);

maxItr=4e3; thresh=1e-10;
u=1;
pOut=[];
prj_C=@(x) max(x,0);

for i=1:5

    a=randn(size(OPT.trueX))*10^(i-3)+OPT.trueX;

    opt.outLevel=1;
    opt.debugLevel=2;
    opt=[];
    opt.Lip=@(u)u^2;
    opt.initStep='fixed';
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
    opt.outLevel=1;
    opt=[];
    tStart=tic;
    sp3=sparseProximal(Psi,Psit,prj_C,'admm',opt);
    for j=1:1
        init=[];
        [x3,itr3,pOut3]=sp3.prox(a,u,thresh,maxItr,init);
    end
    t3=toc(tStart);


    fprintf('objective=%10.10g\n', 0.5*sqrNorm(x1-a)+u*sp3.val(x1));
    fprintf('objective=%10.10g\n', 0.5*sqrNorm(x2-a)+u*sp3.val(x2));
    fprintf('objective=%10.10g\n', 0.5*sqrNorm(x3-a)+u*sp3.val(x3));

    fprintf('results: t1=%g, t2=%g, t3=%g\n', t1, t2, t3);

    keyboard
end

