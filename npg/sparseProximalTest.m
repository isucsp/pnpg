clear all; close all;
if(~exist('seed','var')) seed=0; end
RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',seed));

uArray={100,100,1000, 1e3};

maxItr=2e3; thresh=1e-5;
tvTypeArray={'l1'; 'iso'};
aArray={rand(1,120)', (1:12), rand(100,100), rand(2e2,2e2)};

for i=2:length(tvTypeArray)
    for j=1:length(aArray)

        tvType=tvTypeArray{i};
        a=aArray{j};
        u=uArray{j};

        tv1=tvProximal(tvType);
        tStart=tic;
        [x1,itr1,pOut1,out1]=tv1.op(a,u,thresh,maxItr,[]);
        t1=toc(tStart);

        continue

        tStart=tic;
        pars.print = false;
        pars.tv =tvType;
        pars.MAXITER = maxItr;
        pars.epsilon = thresh; 
        pars.init=zeros(size(a));
        [x2,itr2,out2]=denoise_bound_mod(a,u,-inf,inf,pars);
        t2=toc(tStart);

        opt.debugLevel=1;
        opt.debugLevel=0;
        opt.adaptiveStep=false;
        opt=[];
        opt.outLevel=1;
        opt.debugLevel=2;
        tv3=sparseProximal(tvType,[],[],opt);
        tStart=tic;
        [x3,itr3,pOut3,out3]=tv3.op(a,u,thresh,maxItr,[]);
        t3=toc(tStart);


        fprintf('objective=%g\n', 0.5*sqrNorm(x1-a)+u*tv3.penalty(x1));
        fprintf('objective=%g\n', 0.5*sqrNorm(x2-a)+u*tv3.penalty(x2));
        fprintf('objective=%g\n', 0.5*sqrNorm(x3-a)+u*tv3.penalty(x3));

        fprintf('results: t1=%g, t2=%g, t3=%g\n', t1, t2, t3);

        pause;
    end
end


