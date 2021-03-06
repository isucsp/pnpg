clear all;
if(~exist('seed','var')) seed=0; end
RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',seed));

uArray={1,1,10, 1e1};

maxItr=4e3; thresh=1e-10;
tvTypeArray={'l1'; 'iso'};
aArray={rand(1,120)', (1:12), rand(100,100), rand(2e2,2e2)};

for i=2:length(tvTypeArray)
    for j=1:length(aArray)

        tvType=tvTypeArray{i};
        a=aArray{j};
        u=uArray{j};

        opt.debugLevel=2;
        opt.adaptiveStep=false;
        opt.outLevel=1;
        opt=[];
        tv1=tvProximal(tvType,[],'pnpg',opt);
        tStart=tic;
        [x1,itr1,pOut1,out1]=tv1.prox(a,u,thresh,maxItr,[]);
        t1=toc(tStart); opt=[];

        tStart=tic;
        opt.debugLevel=2;
        opt.outLevel=1;
        opt=[];
        pars.print = out1.opt.debugLevel>0;
        pars.tv =tvType;
        pars.MAXITER = maxItr;
        pars.epsilon = thresh; 
        pars.init=zeros(size(a));
        %[x2,itr2,pOut2]=denoise_bound_mod(a,u,-inf,inf,pars);
        %out2=out1;
        tv2=tvProximal(tvType,[],'beck',opt);
        [x2,itr2,pOut2,out2]=tv2.prox(a,u,thresh,maxItr,[]);
        t2=toc(tStart); opt=[];

        opt.debugLevel=1;
        opt.debugLevel=0;
        opt.debugLevel=2;
        opt.outLevel=1;
        opt=[];
        tv3=tvProximal(tvType,[],'npg',opt);
        tStart=tic;
        [x3,itr3,pOut3,out3]=tv3.prox(a,u,thresh,maxItr,[]);
        t3=toc(tStart); opt=[];

        f=@(x)0.5*sqrNorm(x-a)+u*tv3.val(x);
        fprintf('objective=%10.8g, itr=%d, gap=%10.8g\n',f(x1),itr1,out1.gap);
        fprintf('objective=%10.8g, itr=%d, gap=%10.8g\n',f(x2),itr2,out2.gap);
        fprintf('objective=%10.8g, itr=%d, gap=%10.8g\n',f(x3),itr3,out3.gap);

        fprintf('results: t1=%g, t2=%g, t3=%g\n\n', t1, t2, t3);
    end
end

