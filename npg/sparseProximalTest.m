clear all; close all;

uArray={100,100,1000};

maxItr=1e3; thresh=1e-13;
tvTypeArray={'l1'; 'iso'};
aArray={(1:12)', (1:12), rand(100,100)};

for i=length(tvTypeArray)
    for j=1:length(aArray)

        tvType=tvTypeArray{i};
        a=aArray{j};
        u=uArray{j};

        tic;
        pars.print = false;
        pars.tv =tvType;
        pars.MAXITER = maxItr;
        pars.epsilon = thresh; 
        pars.init=zeros(size(a));
        [D,iter,fun_all]=denoise_bound_mod(a,u,-inf,inf,pars);
        t1=toc;
        fprintf('objective=%g\n', 0.5*sqrNorm(D-a)+u*tlv(D,tvType));

        tic;
        tv=sparseProximal(tvType);
        [dtv,itr]=tv.op(a,u,thresh,maxItr);
        t2=toc;
        fprintf('objective=%g\n', 0.5*sqrNorm(dtv-a)+u*tlv(dtv,tvType));

        fprintf('results: diff=%g, t1=%g, t2=%g\n', norm(D-dtv,'fro'), t1, t2);

        pause;
    end
end


