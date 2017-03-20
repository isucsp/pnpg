function PET_TV_Ex(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%
%
%                  PET example, with background noise b
%      Vary the total counts of the measurements, with continuation

if(~exist('op','var')) op='run'; end

switch lower(op)
case 'run'
    filename = [mfilename '.mat'];
    if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
    clear('OPT','C','proximal');
    filename = [mfilename '.mat'];
    OPT.mask=[]; OPT.outLevel=1;
    OPT.maxItr=1e3; OPT.thresh=1e-6; OPT.debugLevel=2; OPT.noiseType='poisson';
    OPT.maxItr=1e4; OPT.thresh=1e-9; OPT.debugLevel=2; OPT.noiseType='poisson';
    C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
    tvType='l1';
    tvType='iso';
    proximal=tvProximal(tvType,C.prox,'pnpg');
    opt=[]; opt.dualGap=true;
    proximal_dualInnerCrit=tvProximal(tvType, C.prox,'pnpg',opt);

    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    K=1;

    as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
    a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];
    atv= [-0.5,-0.5,  0,   0, 0.5,   1];

    for k=1:K
    for i=4:6
    %for i=[1,2,3,5,6]
    %for i = 6;
        j=1;
        [y,Phi,Phit,~,~,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
        NLL=@(x) Utils.poissonModel(x,Phi,Phit,y,OPT.bb);

        fbp{i,1,k}.x=fbpfunc(y);
        fbp{i,1,k}.RMSE=sqrNorm(maskFunc(fbp{i,1,k}.x,OPT.mask)-OPT.trueX)/sqrNorm(OPT.trueX);

        fprintf('fbp RMSE=%f\n',fbp{i,1,k}.RMSE);
        fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
        u_max=1;
        OPT.u = 10^atv(i)*u_max; OPT.proximal=['tv' tvType];

        initSig=C.prox(maskFunc(fbp{i,1,k}.x,OPT.mask));

        fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);

        opt=OPT;
        opt.grad1 = @(y)[diff(y,1,2), zeros(size(y,1),1)];
        opt.grad2 = @(y)[diff(y); zeros(1,size(y,2))];
        opt.div   = @(x1,x2)([-x1(:,1), -diff(x1(:,1:end-1),1,2), x1(:,end-1)] + [-x2(1,:);-diff(x2(1:end-1,:)); x2(end-1,:)]);
        opt.P = 1e10; opt.p = 2;
        opt.alpha_min = 1e-5; opt.alpha_max = 1e2;
        opt.inn_ini  = 1;
        opt.eta = 1e-6;
        opt.thresh=opt.thresh/1e2;
        vmila{i,j,k} = VMILA(y, Phi, Phit, opt.bb,...
            opt.u, opt.grad1, opt.grad2, opt.div, opt.maxItr,...
            initSig, opt.debugLevel>0, {opt.trueX}, opt.eta, opt.P, opt.p,...
            opt.alpha_min, opt.alpha_max, opt.inn_ini,opt.thresh);

        mysave
        continue;
        opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-5;
        tfocs_200_m5 {i,j,k}=Wrapper.tfocs    (Phi,Phit,[],[],y,initSig,opt);

        opt=OPT;
        spiralTV{i,j,k}=Wrapper.SPIRAL (Phi,Phit,[],[],y,initSig,opt);
        mysave

        opt=OPT; opt.thresh=opt.thresh/100;   opt.maxItr=3e3; opt.xxx=pnpg_{i}.cost(end);
        sigma =[ 1e1,1e-2,1e-3, 10^-4,1e-5,1e-7];
        sigma1=[   1,   1,   1, 10^-0,   1,1e-0];
        tau=   [1e-1,1e-1,1e-2, 10^-2,1e-3,1e-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv2 {i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);
        mysave;
        continue

        opt=OPT;
        pnpg_ {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        mysave;


        opt=OPT; opt.thresh=1e-10;
        spiralTV_Long=Wrapper.SPIRAL (Phi,Phit,[],[],y,initSig,opt);

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_d{i,j}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1; opt.epsilonDecRate=0.6;
        pnpg_d_06{i,j}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        OPT.stepIncre=0.5; OPT.stepShrnk=0.5;
        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_d55{i,j}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1; opt.epsilonDecRate=0.6;
        pnpg_d55_06{i,j}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        opt=OPT;
        pnpg_55 {i,j}=pnpg(NLL,proximal,initSig,opt);
        mysave;
        return;

        OPT.stepIncre=0.5; OPT.stepShrnk=0.5;
        opt=OPT;
        pnpg_55   {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0;
        pnpgG5A055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4;
        pnpgG5Aq55{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0;
        pnpgGfA055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4;
        pnpgGfAq55{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.b=0;
        pnpgA055  {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
        pnpg_n0m055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.adaptiveStep=false;
        pnpg_nInf55{i,j,k}=pnpg(NLL,proximal,initSig,opt);

        mysave;

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_d55   {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgG5A0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgG5Aq_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgGfA0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgGfAq_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgA0_d55  {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_n0m0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.adaptiveStep=false; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_nInf_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        mysave

        OPT.stepIncre=0.8; OPT.stepShrnk=0.8;
        opt=OPT;
        pnpg_88   {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0;
        pnpgG5A088{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4;
        pnpgG5Aq88{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0;
        pnpgGfA088{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4;
        pnpgGfAq88{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.b=0;
        pnpgA088  {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
        pnpg_n0m088{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.adaptiveStep=false;
        pnpg_nInf88{i,j,k}=pnpg(NLL,proximal,initSig,opt);

        mysave;

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_d88   {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgG5A0_d88{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgG5Aq_d88{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgGfA0_d88{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgGfAq_d88{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgA0_d88  {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_n0m0_d88{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.adaptiveStep=false; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_nInf_d88{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        mysave

        return


        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        opt.stepShrnk=0.8; opt.stepIncre=0.8;
        pnpg_d88   {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        mysave
        continue;

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        opt.stepShrnk=0.5; opt.stepIncre=0.5;
        pnpg_d55   {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        mysave
        continue;

        opt=OPT; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_d   {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgG5A0_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgG5Aq_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgGfA0_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgGfAq_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.b=0; opt.dualGap=true; opt.relInnerThresh=1;
        pnpgA0_d  {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_n0m0_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.adaptiveStep=false; opt.dualGap=true; opt.relInnerThresh=1;
        pnpg_nInf_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        mysave
        continue;


        opt=OPT; opt.thresh=opt.thresh/100;  % opt.maxItr=1.5e3; opt.xxx=pnpg_{i}.cost(end);
        sigma=[0,0,0,1e-3,10^-4];
        tau=[1,0.9,0.9,1e-3,10^-4];
        opt.sigma=sigma(i); opt.tau=tau(i);
        cptv1 {i,j,k}=CP_TV(Phi,Phit,y,1,tvType,C,initSig,opt);

        continue;


        mysave

        if(i==5) continue; end

        opt=OPT;
        pnpg_   {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0;
        pnpgG5A0{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4;
        pnpgG5Aq{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0;
        pnpgGfA0{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4;
        pnpgGfAq{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.b=0;
        pnpgA0  {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
        pnpg_n0m0{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.adaptiveStep=false;
        pnpg_nInf{i,j,k}=pnpg(NLL,proximal,initSig,opt);

        mysave;

%       % following are methods for weighted versions
%       ty=max(sqrt(y),1);
%       wPhi=@(xx) Phi(xx)./ty;
%       wPhit=@(xx) Phit(xx./ty);
%       wy=(y-opt.bb(:))./ty;
%       wu_max=pNorm([](wPhit(wy)),inf);
%       opt.noiseType='gaussian';
    end
    end

case lower('tspAddition')
    filename = [mfilename '.mat']; load(filename);
    fprintf('PET Poisson TV example for TSP\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    testList={
        'pnpg_d',
        'pnpg_d_06',
        'pnpg_',
        'pnpg_d55',
        'pnpg_d55_06',
        'pnpg_d55'};

    pnpgd55List={
        'pnpg_d55',
        'pnpgG5A0_d55',
        'pnpgG5Aq_d55',
        'pnpgGfA0_d55',
        'pnpgGfAq_d55',
        'pnpgA0_d55',
        'pnpg_n0m0_d55',
        'pnpg_nInf_d55'};

    pnpgd88List={
        'pnpg_d88',
        'pnpgG5A0_d88',
        'pnpgG5Aq_d88',
        'pnpgGfA0_d88',
        'pnpgGfAq_d88',
        'pnpgA0_d88',
        'pnpg_n0m0_d88',
        'pnpg_nInf_d88'};

    pnpg55List={
        'pnpg_55',
        'pnpgG5A055',
        'pnpgG5Aq55',
        'pnpgGfA055',
        'pnpgGfAq55',
        'pnpgA055',
        'pnpg_n0m055',
        'pnpg_nInf55'};

    pnpg88List={
        'pnpg_88',
        'pnpgG5A088',
        'pnpgG5Aq88',
        'pnpgGfA088',
        'pnpgGfAq88',
        'pnpgA088',
        'pnpg_n0m088',
        'pnpg_nInf88'};

    pnpgdList={
        'pnpg_d'};

    pnpgList={
        'pnpg_',
        'pnpgG5A0',
        'pnpgG5Aq',
        'pnpgGfA0',
        'pnpgGfAq',
        'pnpgA0',
        'pnpg_n0m0',
        'pnpg_nInf'};

    tfocsList={
        'tfocs_200_m6 ',
        'tfocs_200_m9 ',
        'tfocs_200_m12'};

    i=4;
    mc=+inf;
    [mc,pnpgd55Var]=minAndName(pnpgd55List,i,mc);
    [mc,pnpgd88Var]=minAndName(pnpgd88List,i,mc);
    [mc,pnpg55Var ]=minAndName(pnpg55List ,i,mc);
    [mc,pnpg88Var ]=minAndName(pnpg88List ,i,mc);
    [mc,pnpgdVar  ]=minAndName(pnpgdList  ,i,mc);
    [mc,pnpgVar   ]=minAndName(pnpgList   ,i,mc);
    [mc,testVar   ]=minAndName(testList   ,i,mc);

    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgd55Var{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgd88Var{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpg55Var{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpg88Var{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgdVar{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgVar{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),testVar{:});


    spiral=spiralTV;
    tfocs=tfocs_200_m5;

    % mIdx=6 is also good
    mIdx=4; as=1; k=1;
    mc=min([  min(    pnpg_{mIdx,as,k}.cost)
              min(  pnpg_n0{mIdx,as,k}.cost)
              min(   spiral{mIdx,as,k}.cost)
              min(pnpg_nInf{mIdx,as,k}.cost)
              min(pnpg_n0m0{mIdx,as,k}.cost)
              min(    tfocs{mIdx,as,k}.cost)
              min(   pnpgA0{mIdx,as,k}.cost)
              min( pnpgG5A0{mIdx,as,k}.cost)
              min( pnpgG5Aq{mIdx,as,k}.cost)
              min( pnpgGfA0{mIdx,as,k}.cost)
              min( pnpgGfAq{mIdx,as,k}.cost)
              min(    vmila{mIdx,as,k}.cost)
              min(    cptv1{mIdx,as,k}.cost)
              min(    cptv2{mIdx,as,k}.cost) ]);

    fields_={'RMSE','time','cost'};
    forSave=addTrace(       pnpg_{mIdx,as,k},     [],fields_,mc); %  1- 4
    forSave=addTrace(     pnpg_n0{mIdx,as,k},forSave,fields_,mc); %  5- 8
    forSave=addTrace(      spiral{mIdx,as,k},forSave,fields_,mc); %  9-12
    forSave=addTrace(   pnpg_nInf{mIdx,as,k},forSave,fields_,mc); % 13-16
    forSave=addTrace(   pnpg_n0m0{mIdx,as,k},forSave,fields_,mc); % 17-20
    forSave=addTrace(       tfocs{mIdx,as,k},forSave,fields_,mc); % 21-24
    forSave=addTrace(      pnpgA0{mIdx,as,k},forSave,fields_,mc); % 25-28
    forSave=addTrace(    pnpgG5A0{mIdx,as,k},forSave,fields_,mc); % 29-32
    forSave=addTrace(    pnpgG5Aq{mIdx,as,k},forSave,fields_,mc); % 33-26
    forSave=addTrace(    pnpgGfA0{mIdx,as,k},forSave,fields_,mc); % 37-40
    forSave=addTrace(    pnpgGfAq{mIdx,as,k},forSave,fields_,mc); % 41-44
    save('cost_itrPET_TV.data','forSave','-ascii');

    o=vmila{mIdx,as,k};
    vmila{mIdx,as,k}.RMSE=o.RMSE{1}(1:length(o.cost));
    fields_={'RMSE','time','cost'};
    forSave=addTrace(       vmila{mIdx,as,k},     [],fields_,mc); %  1- 4
    forSave=addTrace(       cptv1{mIdx,as,k},forSave,fields_,mc); %  5- 8
    forSave=addTrace(       cptv2{mIdx,as,k},forSave,fields_,mc); %  9-12
    forSave=addTrace(      pnpg_d{mIdx,as,k},forSave,fields_,mc); % 12-16
    save('cost_itrPET_TV_1.data','forSave','-ascii');
    paperDir='~/research/myPaper/asilomar2014/';
    system(['mv cost_itrPET_TV.data cost_itrPET_TV_1.data ' paperDir]);

    nn=128;
    xtrue = read_zubal_emis('nx', nn, 'ny', nn);
    idx=5;
    fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf(' tfocs: %g%%\n', tfocs{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
    img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
    img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_TV_pet.png')
    img=showImgMask( tfocs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'tfocs_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'tfocs_TV_pet.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_TV_pet.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_TV_pet.png')

    idx=4;
    fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf(' tfocs: %g%%\n', tfocs{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
    img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
    img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_TV4_pet.png')
    img=showImgMask( tfocs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'tfocs_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'tfocs_TV4_pet.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_TV4_pet.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_TV4_pet.png')

case lower('plotTV')
    filename = [mfilename '.mat']; load(filename);
    fprintf('PET Poisson TV example\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    K = 1:1;

    npgTVcTime= mean(Cell.getField(npgTVc(:,1,K),'time'),3);
    npgTVcCost= mean(Cell.getField(npgTVc(:,1,K),'cost'),3);
    npgTVcRMSE= mean(Cell.getField(npgTVc(:,1,K),'RMSE'),3);

    npgTVTime = mean(Cell.getField(npgTV (:,1,K),'time'),3);
    npgTVCost = mean(Cell.getField(npgTV (:,1,K),'cost'),3);
    npgTVRMSE = mean(Cell.getField(npgTV (:,1,K),'RMSE'),3);

    spiralTVTime = mean(Cell.getField(spiralTV (:,1,K),'time'),3);
    spiralTVCost = mean(Cell.getField(spiralTV (:,1,K),'cost'),3);
    spiralTVRMSE = mean(Cell.getField(spiralTV (:,1,K),'RMSE'),3);

    npgsTime  = mean(Cell.getField(  npgs(:,1,K),'time'),3);
    npgsCost  = mean(Cell.getField(  npgs(:,1,K),'cost'),3);
    npgsRMSE  = mean(Cell.getField(  npgs(:,1,K),'RMSE'),3);

    npgscTime = mean(Cell.getField( npgsc(:,1,K),'time'),3);
    npgscCost = mean(Cell.getField( npgsc(:,1,K),'cost'),3);
    npgscRMSE = mean(Cell.getField( npgsc(:,1,K),'RMSE'),3);

    vmilaTime = mean(Cell.getField( vmila(:,1,K),'time'),3);
    vmilaCost = mean(Cell.getField( vmila(:,1,K),'cost'),3);
    vmilaRMSE = mean(Cell.getField( vmila(:,1,K),'RMSE'),3);

    cptv1Time = mean(Cell.getField( cptv1(:,1,K),'time'),3);
    cptv1Cost = mean(Cell.getField( cptv1(:,1,K),'cost'),3);
    cptv1RMSE = mean(Cell.getField( cptv1(:,1,K),'RMSE'),3);

    cptv2Time = mean(Cell.getField( cptv2(:,1,K),'time'),3);
    cptv2Cost = mean(Cell.getField( cptv2(:,1,K),'cost'),3);
    cptv2RMSE = mean(Cell.getField( cptv2(:,1,K),'RMSE'),3);

    fbpRMSE   = mean(Cell.getField(   fbp(:,1,K),'RMSE'),3);

    figure;
    loglog(count,npgTVRMSE,'r-*'); hold on;
    loglog(count,   fbpRMSE,'b-o');
    loglog(count,spiralTVRMSE,'k-^');
    loglog(count,  npgTVcRMSE,'k*-.');
    loglog(count,  npgsRMSE,'c>-');
    loglog(count, npgscRMSE,'gs-');
    legend('npgTV','fbp','spiralTV','npgTVc','npgs','npgsc');

    figure;
    loglog(count,   npgTVTime,'r-*'); hold on;
    loglog(count,spiralTVTime,'k-^');
    loglog(count,  npgTVcTime,'k*-.');
    loglog(count,  npgsTime,'c>-');
    loglog(count, npgscTime,'gs-');
    legend('npgTV','spiralTV','npgTVc','npgs','npgsc');

    forSave=[npgTVTime, npgTVcTime, npgsTime, npgscTime, spiralTVTime,...
        npgTVCost, npgTVcCost, npgsCost, npgscCost, spiralTVCost,...
        npgTVRMSE, npgTVcRMSE, npgsRMSE, npgscRMSE, spiralTVRMSE,...
        fbpRMSE, count(:)];
    save('varyCntPETTV.data','forSave','-ascii');

    forSave=[]; t=0; mIdx=5; k=1;
    out=   npgTV_n4;
    t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=spiralTV_Long;
    t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=pnpgTV_noAdpStpLong;
    t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=   pnpgTV_n1;
    t=t+1; forSave(1:length(out.stepSize),t)=out.stepSize;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    save('stepSize.data','forSave','-ascii');

    figure; semilogy(forSave(:,1),'r'); hold on;
    semilogy(forSave(:,3),'g');
    semilogy(forSave(:,5),'b');
    title('step size versus number of iterations');
    legend('npgTV','spiralTV','npgTV noAdaptive Step');

    forSave=[]; t=0; mIdx=5; k=1;
    out=   npgTV_n4;
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=  npgTV_n1;
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=  npgTV_noAdpStpLong;
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    t=t+1; forSave(1:length(out.difX),t)=out.difX;
    out=  spiralTV_Long;
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    t=t+1; forSave(1:length(out.difX),t)=out.difX;

    save('cost_itrPETTV.data','forSave','-ascii');
    mincost=reshape(forSave(:,[1,4,7,11]),[],1); 
    mincost=min(mincost(mincost~=0));

    figure;
    semilogy(forSave(:,3),forSave(:,1)-mincost,'r'); hold on;
    semilogy(forSave(:,6),forSave(:,4)-mincost,'g');
    semilogy(forSave(:,9),forSave(:,7)-mincost,'b');
    semilogy(forSave(:,13),forSave(:,11)-mincost,'k');
    legend('npgTV n4','npgTV n1','npgTV nInf','spiralTV');
    hold on;
    idx=min(find(forSave(:,10)<1e-6));
    plot(forSave(idx,9),forSave(idx,7)-mincost,'bo');
    xxx=idx;
    idx=min(find(forSave(10:end,14)<1e-6))+10;
    plot(forSave(idx,13),forSave(idx,11)-mincost,'k*');
    xxx=[xxx;idx];  xxx=xxx(:)';
    save('cost_itrPETTVidx.data','xxx','-ascii');

    figure; semilogy(forSave(:,3),forSave(:,2),'r'); hold on;
    semilogy(forSave(:,6),forSave(:,5),'g');
    semilogy(forSave(:,9),forSave(:,8),'b');
    semilogy(forSave(:,13),forSave(:,12),'k');
    legend('npgTV n4','npgTV n1','npgTV nInf','spiralTV');

    keyboard

    nn=128;
    xtrue = read_zubal_emis('nx', nn, 'ny', nn);
    % attenuation map
    mumap = read_zubal_attn('nx', nn, 'ny', nn);
    imwrite(xtrue/max(xtrue(:)),'pet.png');
    imwrite(mumap/max(mumap(:)),'mumap.png');

    idx=5;
    fprintf('   NPGTV: %g%%\n',   npgTV{idx}.RMSE(end)*100);
    fprintf('SPIRALTV: %g%%\n',spiralTV{idx}.RMSE(end)*100);
    fprintf('     FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
    fprintf('    NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
    img=npg{idx}.x; mask=npg{idx}.opt.mask;
    img=showImgMask(   npgTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPGTV_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPGTV_pet.png')
    img=showImgMask(spiralTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRALTV_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRALTV_pet.png')
    img=showImgMask(     fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,     'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),     'FBP_pet.png')
    img=showImgMask(    npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,    'NPGs_pet.eps','psc2'); imwrite(img/max(xtrue(:)),    'NPGs_pet.png')

    idx=4;
    fprintf('   NPGTV: %g%%\n',   npgTV{idx}.RMSE(end)*100);
    fprintf('SPIRALTV: %g%%\n',spiralTV{idx}.RMSE(end)*100);
    fprintf('     FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
    fprintf('    NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
    img=npg{idx}.x; mask=npg{idx}.opt.mask;
    img=showImgMask(   npgTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPGTV_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPGTV_pet2.png')
    img=showImgMask(spiralTV{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRALTV_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRALTV_pet2.png')
    img=showImgMask(     fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,     'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),     'FBP_pet2.png')
    img=showImgMask(    npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,    'NPGs_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),    'NPGs_pet2.png')

    decide=input(sprintf('start to copy to %s [y/N]?',paperDir));
    if strcmpi(decide,'y')
        system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
    end
    system('rm *_pet.eps *_pet2.eps *_pet2.png');
    close all;

case 'plot'
    filename = [mfilename '.mat'];
    load(filename);
    fprintf('PET Poisson example\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    K = 1:1;

    npgTime      = mean(Cell.getField(npg     (:,1,K),'time'),3);
    npgCost      = mean(Cell.getField(npg     (:,1,K),'cost'),3);
    npgRMSE      = mean(Cell.getField(npg     (:,1,K),'RMSE'),3);
    npgcTime     = mean(Cell.getField(npgc    (:,1,K),'time'),3);
    npgcCost     = mean(Cell.getField(npgc    (:,1,K),'cost'),3);
    npgcRMSE     = mean(Cell.getField(npgc    (:,1,K),'RMSE'),3);
    npgsTime     = mean(Cell.getField(npgs    (:,1,K),'time'),3);
    npgsCost     = mean(Cell.getField(npgs    (:,1,K),'cost'),3);
    npgsRMSE     = mean(Cell.getField(npgs    (:,1,K),'RMSE'),3);
    npgscTime    = mean(Cell.getField(npgsc   (:,1,K),'time'),3);
    npgscCost    = mean(Cell.getField(npgsc   (:,1,K),'cost'),3);
    npgscRMSE    = mean(Cell.getField(npgsc   (:,1,K),'RMSE'),3);
    spiralTime   = mean(Cell.getField(spiral  (:,1,K),'time'),3);
    spiralCost   = mean(Cell.getField(spiral  (:,1,K),'cost'),3);
    spiralRMSE   = mean(Cell.getField(spiral  (:,1,K),'RMSE'),3);

    fbpRMSE      = mean(Cell.getField(fbp     (:,1,K),'RMSE'),3);

    figure;
    loglog(count,     npgRMSE,'r-*'); hold on;
    loglog(count,     fbpRMSE,'b-o');
    loglog(count,  spiralRMSE,'k-^');
    loglog(count,    npgcRMSE,'k*-.');
    loglog(count,    npgsRMSE,'c>-');
    loglog(count,   npgscRMSE,'gs-');
    legend('npg','fbp','spiral','npgc','npgs','npgsc');

    figure;
    loglog(count,     npgTime,'r-*'); hold on;
    loglog(count,  spiralTime,'k-^');
    loglog(count ,   npgcTime,'k*-.');
    loglog(count,    npgsTime,'c>-');
    loglog(count,   npgscTime,'gs-');
    legend('npg','spiral','npgc','npgs','npgsc');

    forSave=[npgTime, npgcTime, npgsTime, npgscTime, spiralTime,...
        npgCost, npgcCost, npgsCost, npgscCost, spiralCost,...
        npgRMSE, npgcRMSE, npgsRMSE, npgscRMSE, spiralRMSE,...
        fbpRMSE, count(:)...
        ];
    save('varyCntPET.data','forSave','-ascii');

    forSave=[]; t=0; mIdx=5; k=1;
    out=  npgc{mIdx,1,k};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=   npg{mIdx,1,k};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=spiral{mIdx,1,k};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=npgs{mIdx,1,k};
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;
    out=npgsc_s;
    t=t+1; forSave(1:length(out.cost),t)=out.cost;
    t=t+1; forSave(1:length(out.RMSE),t)=out.RMSE;
    t=t+1; forSave(1:length(out.time),t)=out.time;

    save('cost_itrPET.data','forSave','-ascii');
    mincost=reshape(forSave(:,[1,4,7]),[],1); 
    mincost=min(mincost(mincost~=0));

    figure;
    semilogy(forSave(:,3),forSave(:,1)-mincost,'r'); hold on;
    semilogy(forSave(:,6),forSave(:,4)-mincost,'g');
    semilogy(forSave(:,9),forSave(:,7)-mincost,'b');
    if(mIdx==5 && k==1)
        semilogy(forSave(:,15),forSave(:,13)-min(max(forSave(:,13),0)),'c:');
        legend('npgc','npg','spiral','npgsc');
    else
        legend('npgc','npg','spiral');
    end
    figure; semilogy(forSave(:,3),forSave(:,2),'r'); hold on;
    semilogy(forSave(:,6),forSave(:,5),'g'); semilogy(forSave(:,9),forSave(:,8),'b');
    legend('npgc','npg','spiral');

    keyboard

    nn=128;
    xtrue = read_zubal_emis('nx', nn, 'ny', nn);
    % attenuation map
    mumap = read_zubal_attn('nx', nn, 'ny', nn);
    imwrite(xtrue/max(xtrue(:)),'pet.png');
    imwrite(mumap/max(mumap(:)),'mumap.png');

    idx=5;
    fprintf('   NPG: %g%%\n',   npg{idx}.RMSE(end)*100);
    fprintf('  NPGc: %g%%\n',  npgc{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
    fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
    fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
    fprintf('NPGscS: (%g%%, %g%%)\n',    npgsc_s.RMSE(end)*100,rmseTruncate(npgsc_s   )*100);
    img=npg{idx}.x; mask=npg{idx}.opt.mask;
    img=showImgMask(   npg{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPG_pet.png')
    img=showImgMask(  npgc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGc_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGc_pet.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet.png')
    img=showImgMask(  npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGs_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGs_pet.png')
    img=showImgMask( npgsc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGsc_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGsc_pet.png')
    img=showImgMask( npgsc_s.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGscS_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGscS_pet.png')

    idx=4;
    fprintf('   NPG: %g%%\n',   npg{idx}.RMSE(end)*100);
    fprintf('  NPGc: %g%%\n',  npgc{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},npg{idx}.opt.trueX)*100);
    fprintf('  NPGs: (%g%%, %g%%)\n',  npgs{idx}.RMSE(end)*100,rmseTruncate( npgs{idx})*100);
    fprintf(' NPGsc: (%g%%, %g%%)\n', npgsc{idx}.RMSE(end)*100,rmseTruncate(npgsc{idx})*100);
    img=npg{idx}.x; mask=npg{idx}.opt.mask;
    img=showImgMask(   npg{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'NPG_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'NPG_pet2.png')
    img=showImgMask(  npgc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGc_pet2.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet2.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet2.png')
    img=showImgMask(  npgs{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'NPGs_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'NPGs_pet2.png')
    img=showImgMask( npgsc{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'NPGsc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)), 'NPGsc_pet2.png')

    keyboard
    paperDir='~/research/myPaper/asilomar2014/'
    system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
    system('rm *_pet.eps *_pet2.eps *_pet2.png');
    close all;
case 'fullplot'
    filename = [mfilename '.mat'];
    load(filename);

    k=1;
    aa =(3:-0.5:-6);
    for i=1:length(count)
        npgContRMSE  {i,k} = [  npgFull{i,k}.contRMSE(:);  npgFull{i,k}.RMSE(end)]; out=npgContRMSE{i,k};
        fprintf('i=%d, good a = 1e%g NPG\n',i,max((aa(out==min(out)))));
        npgsContRMSE {i,k} = [ npgsFull{i,k}.contRMSE(:); npgsFull{i,k}.RMSE(end)]; out=npgsContRMSE{i,k};
        fprintf('i=%d, good a = 1e%g NPGs\n',i,max((aa(out==min(out)))));
        npgTVContRMSE{i,k} = [npgTVFull{i,k}.contRMSE(:);npgTVFull{i,k}.RMSE(end)]; out=npgTVContRMSE{i,k};
        fprintf('i=%d, good a = 1e%g NPG_TV\n',i,max((aa(out==min(out)))));
    end

    for i=1:length(count)
        figure;
        semilogy(aa(1:length(npgContRMSE{i})),npgContRMSE{i},'r-*'); hold on;
        semilogy(aa(1:length(npgsContRMSE{i})),npgsContRMSE{i},'g-o');
        semilogy(aa(1:length(npgTVContRMSE{i})),npgTVContRMSE{i},'b-s');
        title(num2str(i));
        legend('NPG','NPGs','NPG-TV');
        aaa(i)=min(npgContRMSE{i});
        bbb(i)=min(npgsContRMSE{i});
        ccc(i)=min(npgTVContRMSE{i});
    end
    figure; semilogy(aaa,'r-*'); hold on;
    semilogy(bbb,'g-o');
    semilogy(ccc,'b-s');
    title('rmse vs count');
    legend('NPG','NPGs','NPG-TV');


    keyboard

    K=1;
    m = [1e4 1e5 1e6 1e7 1e8 1e9];
    fprintf('Poisson example\n');

    npgTime    = mean(Cell.getField(    npg(:,:,1:K),'time'),3);
    npgcTime   = mean(Cell.getField(   npgc(:,:,1:K),'time'),3);
    npgsTime   = mean(Cell.getField(   npgs(:,:,1:K),'time'),3);
    npgscTime  = mean(Cell.getField(  npgsc(:,:,1:K),'time'),3);
    spiralTime = mean(Cell.getField( spiral(:,:,1:K),'time'),3);

    npgCost    = mean(Cell.getField(    npg(:,:,1:K),'cost'),3);
    npgcCost   = mean(Cell.getField(   npgc(:,:,1:K),'cost'),3);
    npgsCost   = mean(Cell.getField(   npgs(:,:,1:K),'cost'),3);
    npgscCost  = mean(Cell.getField(  npgsc(:,:,1:K),'cost'),3);
    spiralCost = mean(Cell.getField( spiral(:,:,1:K),'cost'),3);

    npgRMSE    = mean(Cell.getField(    npg(:,:,1:K),'RMSE'),3);
    npgcRMSE   = mean(Cell.getField(   npgc(:,:,1:K),'RMSE'),3);
    npgsRMSE   = mean(Cell.getField(   npgs(:,:,1:K),'RMSE'),3);
    npgscRMSE  = mean(Cell.getField(  npgsc(:,:,1:K),'RMSE'),3);
    spiralRMSE = mean(Cell.getField( spiral(:,:,1:K),'RMSE'),3);

    aIdx=4;
    figure;
    loglog(m,    npgRMSE(:,aIdx),'r-*'); hold on;
    loglog(m,   npgsRMSE(:,aIdx),'c-p');
    loglog(m, spiralRMSE(:,aIdx),'k-^');
    loglog(m,   npgcRMSE(:,aIdx),'k*-.');
    loglog(m,  npgscRMSE(:,aIdx),'bs-.');
    legend('npg','npgs','spiral','npgc','npgsc');

    figure;
    loglog(m,    npgTime(:,aIdx),'r-*' ); hold on;
    loglog(m,   npgsTime(:,aIdx),'c-p' );
    loglog(m, spiralTime(:,aIdx),'k-^' );
    loglog(m,   npgcTime(:,aIdx),'k*-.');
    loglog(m,  npgscTime(:,aIdx),'bs-.');
    legend('npg','npgs','spiral','npgc','npgsc');

    forSave=[npgTime(:,aIdx), npgsTime(:,aIdx), npgcTime(:,aIdx), npgscTime(:,aIdx), spiralTime(:,aIdx),...
        npgCost(:,aIdx), npgsCost(:,aIdx), npgcCost(:,aIdx), npgscCost(:,aIdx), spiralCost(:,aIdx),...
        npgRMSE(:,aIdx), npgsRMSE(:,aIdx), npgcRMSE(:,aIdx), npgscRMSE(:,aIdx), spiralRMSE(:,aIdx),...
        m(:)];
    save('varyMeasurementPoisson.data','forSave','-ascii');

    % time cost RMSE
    forSave=[count(:),meanOverK(   fbp,'RMSE'),...
        meanOverK(    pnpg_),...
        meanOverK(  pnpg_n0),...
        meanOverK(   spiral),...
        meanOverK(    tfocs),...
        ];
    save('varyCntPET_TV.data','forSave','-ascii');


end

pnpg_1_3=[];
spiral_m5=[];
spiral_m6=[];
pnpgTV_noAdpStpLong=[];
pnpg_d_06=[];
pnpg_d55_06=[];

function [mc,varList] = minAndName(nameList,i,mc)
    if(~exist('mc','var'))
        mc=+inf;
    end
    for ii=1:length(nameList)
        a=eval(nameList{ii});
        mc=min(mc,min(a{i}.cost(:)));
        a{i}.name=nameList{ii};
        varList{ii}=a{i};
    end
end


end

function [a,b,c]=meanOverK(method,field)
    if(nargin==2)
        a=mean(Cell.getField(method,field),3);
    else
        a=mean(Cell.getField(method,'time'),3);
        b=mean(Cell.getField(method,'cost'),3);
        c=mean(Cell.getField(method,'RMSE'),3);
        a=[a b c];
    end
end
function forSave=addTrace(method,forSave,fields,mc)
    len=1000;
    tt=method.(fields{1});
    itr=linspace(1,length(tt),len);
    if(~exist('fields','var'))
        fields={'time','cost','RMSE'};
    end
    for i=1:length(fields);
        tt=method.(fields{i});
        if(iscell(tt) && length(tt)==1)
            tt=tt{1};
        end
        if(strcmpi(fields{i},'cost') && exist('mc','var'))
            tt=(tt-mc)/mc;
        end
        ss=interp1(1:length(tt),tt,itr);
        data(:,i)=reshape(ss,[],1);
    end
    data=[itr(:) data];
    forSave=appendColumns(data,forSave);
end
function forSave = appendColumns(col,forSave)
    [r,c]=size(forSave);
    forSave(1:size(col,1),c+1:c+size(col,2))=col;
end



