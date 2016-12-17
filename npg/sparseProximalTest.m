clear all; close all;
if(~exist('seed','var')) seed=0; end
RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',seed));

uArray={100,100,1000, 1e3};

maxItr=2e3; thresh=1e-4;
tvTypeArray={'l1'; 'iso'};
aArray={rand(1,120)', (1:12), rand(100,100), rand(5e2,5e2)};

for i=1:length(tvTypeArray)
    for j=1:length(aArray)

        tvType=tvTypeArray{i};
        a=aArray{j};
        u=uArray{j};

        tStart=tic;
        pars.print = true;
        pars.tv =tvType;
        pars.MAXITER = maxItr;
        pars.epsilon = thresh; 
        pars.init=zeros(size(a));
        [D,iter,fun_all]=denoise_bound_mod(a,u,-inf,inf,pars);
        t1=toc(tStart);

        tStart=tic;
        opt=[];
        opt.debugLevel=2;
        tv=sparseProximal(tvType,[],[],opt);
        [dtv,itr,pInit,out]=tv.op(a,u,thresh,maxItr,[]);
        t2=toc(tStart);
        fprintf('objective=%g\n', 0.5*sqrNorm(D-a)+u*tv.penalty(D));

        fprintf('objective=%g\n', 0.5*sqrNorm(dtv-a)+u*tv.penalty(dtv));

        fprintf('results: diff=%g, t1=%g, t2=%g\n', norm(D-dtv,'fro'), t1, t2);

        pause;
    end
end


