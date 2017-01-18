function PET_DWT_Ex(op)
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
    clear('OPT','C','proximal','PROXOPT');
    filename = [mfilename '.mat'];
    OPT.mask=[]; OPT.outLevel=1;
    OPT.maxItr=2e3; OPT.thresh=1e-9; OPT.debugLevel=2; OPT.noiseType='poisson';
    OPT.maxItr=2e4; OPT.thresh=1e-9; OPT.debugLevel=2; OPT.noiseType='poisson';
    C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);
    PROXOPT.Lip=@(u)u^2; PROXOPT.initStep='fixed';
    PROXOPT.adaptiveStep=false; PROXOPT.backtracking=false;
    proximal=sparseProximal(Psi,Psit,C.prox,'pnpg',PROXOPT);
    proxOpt=PROXOPT;  proxOpt.dualGap=true;
    proximal_dualInnerCrit=sparseProximal(Psi,Psit,C.prox,'pnpg',proxOpt);

    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    K=1;

    as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
    a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];
    atv= [-0.5,-0.5,  0,   0, 0.5,   1];

    for k=1:K
    for i=5
        j=1;
        [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
        NLL=@(x) Utils.poissonModel(x,Phi,Phit,y,OPT.bb);

        fbp{i,1,k}.x=fbpfunc(y);
        fbp{i,1,k}.RMSE=sqrNorm(maskFunc(fbp{i,1,k}.x,OPT.mask)-OPT.trueX)/sqrNorm(OPT.trueX);

        fprintf('fbp RMSE=%f\n',fbp{i,1,k}.RMSE);
        fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
        u_max=1;
        OPT.u = u_max*10.^a(i);

        initSig=C.prox(maskFunc(fbp{i,1,k}.x,OPT.mask));

        fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);

        % BEGIN experiment region,  to delete in the end
        % END experiment region,  to delete in the end

        opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-6;
        opt.maxInnerItr=100;
        tfocs_200_m6_100 {i,j,k}=Wrapper.tfocs(Phi,Phit,Psi,Psit,y,initSig,opt);
        opt=OPT; opt.restartEvery=200; opt.innerThresh=1e-9;
        opt.maxInnerItr=100;
        tfocs_200_m9_100 {i,j,k}=Wrapper.tfocs(Phi,Phit,Psi,Psit,y,initSig,opt);
        mysave;
        return;

        opt=OPT; opt.maxInnerItr=1e3;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpg_55 {i,j}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgG5Aq55{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgG5A055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgGfAq55{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgGfA055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.b=0;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgA055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpg_n055{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.adaptiveStep=false;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpg_nInf55{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        mysave

        opt=OPT; opt.dualGap=true;
        pnpg_d{i,j}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0; opt.dualGap=true;
        pnpgG5A0_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4; opt.dualGap=true;
        pnpgG5Aq_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0; opt.dualGap=true;
        pnpgGfA0_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4; opt.dualGap=true;
        pnpgGfAq_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.b=0; opt.dualGap=true;
        pnpgA0_d  {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false; opt.dualGap=true;
        pnpg_n0_d {i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.adaptiveStep=false; opt.dualGap=true;
        pnpg_nInf_d{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        mysave

        opt=OPT; opt.maxInnerItr=1e3;
        pnpg_ {i,j}=pnpg(NLL,proximal,initSig,opt);
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
        pnpg_n0 {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        opt=OPT; opt.adaptiveStep=false;
        pnpg_nInf{i,j,k}=pnpg(NLL,proximal,initSig,opt);
        mysave;
         
        opt=OPT; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpg_d55{i,j}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=0; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgG5A0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=5; opt.b=1/4; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgG5Aq_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=0; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgGfA0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.gamma=15; opt.b=1/4; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgGfAq_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.b=0; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpgA0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.cumuTol=0; opt.incCumuTol=false; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpg_n0_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);
        opt=OPT; opt.adaptiveStep=false; opt.dualGap=true;
        opt.stepIncre=0.5; opt.stepShrnk=0.5;
        pnpg_nInf_d55{i,j,k}=pnpg(NLL,proximal_dualInnerCrit,initSig,opt);

        % this version maybe able to used to achieve the smallest objective
        opt=OPT;
        proxOpt=PROXOPT;
        opt.proximal=sparseProximal(Psi,Psit,C.prox,'admm',proxOpt);
        pnpg_admm   {i,j,k}=pnpg(NLL,opt.proximal,initSig,opt);

        mysave;

        continue;

        opt=OPT; 
        opt.stepShrnk=0.8; opt.stepIncre=0.8;
        pnpg_88 {i,j}=pnpg(NLL,proximal,initSig,opt);

        continue;
         
        opt=OPT; opt.thresh=opt.thresh/100;  % opt.maxItr=1.5e3;
        L=1e5;                               % opt.xxx=pnpg_{i}.cost(end);
        sigma=[0,0,0,0,10^-4];
        tau=[1,0.9,0.9,0.8,10^-3];
        opt.sigma=sigma(i); opt.tau=1/L/opt.sigma*tau(i);
        cpdwt1 {i,j}=CP_DWT(Phi,Phit,y,1,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;  % opt.maxItr=3e3;
        L=1e5;                               % opt.xxx=pnpg_{i}.cost(end);
        sigma=[0,0,0,0,1e-6];
        sigma1=[0,1e-3,1e-2,1e-1,1];
        tau=[1,1,1,1,10^-3];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt2 {i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.maxItr=10*opt.maxItr;
        opt.L=1/mean(pnpg_{i,j,k}.stepSize);
        sigma=[0,0,0,0,1];
        tau=[0,1,1,1,1];
        rho=[1,1,1,1,1];
        opt.sigma=sigma(i); opt.tau=tau(i); opt.rho=rho(i);
        H.exact=true;
        H.val=@(s) opt.u*norm(s(:),1);
        H.proxConj=@(a,v) max(min(a,opt.u),-opt.u);
        condat {i,j,k}=pds(NLL,H,Psi,Psit,C,opt.L,initSig,opt);

        opt=OPT; opt.innerThresh=1e-6;
        spiral_m6   {i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
        opt=OPT; opt.innerThresh=1e-9;
        spiral_m9   {i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
        opt=OPT; opt.innerThresh=1e-12;
        spiral_m12  {i,j,k}=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt);
        mysave

        opt=OPT;   % opt.maxItr=1e3; opt.debugLevel=1;
        opt.grad1 = @(y)Psit(y);
        opt.grad2 = @(y)0;
        opt.div   = @(x1,x2) Psi(x1);
        opt.P = 1e5; opt.p = 2;
        opt.alpha_min = 1e-5; opt.alpha_max = 1e2;
        opt.inn_ini  = 1;
        opt.eta = 1e-6;
        opt.thresh=opt.thresh/1e2;
        vmila{i,j,k} = VMILA(y, Phi, Phit, opt.bb,...
            opt.u, opt.grad1, opt.grad2, opt.div, opt.maxItr,...
            initSig, opt.debugLevel>0, {opt.trueX}, opt.eta, opt.P, opt.p,...
            opt.alpha_min, opt.alpha_max, opt.inn_ini,opt.thresh);
        mysave;
    end
    end

case lower('plot')
    filename = [mfilename '.mat']; load(filename);
    fprintf('PET Poisson l1 example\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    pnpgd55List={
        'pnpg_d55',
        'pnpgG5A0_d55',
        'pnpgG5Aq_d55',
        'pnpgGfA0_d55',
        'pnpgGfAq_d55',
        'pnpgA0_d55',
        'pnpg_n0_d55',
        'pnpg_nInf_d55'};

    pnpg55List={
        'pnpg_55',
        'pnpgG5A055',
        'pnpgG5Aq55',
        'pnpgGfA055',
        'pnpgGfAq55',
        'pnpgA055',
        'pnpg_n055',
        'pnpg_nInf55'};

    pnpgdList={
        'pnpg_d',
        'pnpgG5A0_d',
        'pnpgG5Aq_d',
        'pnpgGfA0_d',
        'pnpgGfAq_d',
        'pnpgA0_d',
        'pnpg_n0_d',
        'pnpg_nInf_d'};

    pnpgList={
        'pnpg_',
        'pnpgG5A0',
        'pnpgG5Aq',
        'pnpgGfA0',
        'pnpgGfAq',
        'pnpgA0',
        'pnpg_n0',
        'pnpg_nInf'};

    tfocsList={
        'tfocs_200_m6 ',
        'tfocs_200_m9 '};

    spiralList={
        'spiral_m6  ',
        'spiral_m9  ',
        'spiral_m12 ',
    };

    testList={
        'pnpg_d',
        'pnpg_',
        'pnpg_d55',
        'pnpg_55'};

    nameList={
        '    pnpg_',
        '   spiral',
        'pnpg_nInf',
        '  pnpg_n0',
        '    tfocs',
        '   pnpgA0',
        ' pnpgG5A0',
        ' pnpgG5Aq',
        ' pnpgGfA0',
        ' pnpgGfAq',
        '   cpdwt1',
        '   cpdwt1',
        '   pnpg_d'};

    K = 1;
          fbp=      fbp(:,:,1:K);
        pnpg_=    pnpg_(:,:,1:K);
       spiral=   spiral_m5(:,:,1:K);
      pnpg_n0=  pnpg_n0(:,:,1:K);
    pnpg_nInf=pnpg_nInf(:,:,1:K);
        tfocs=tfocs_200_m5(:,:,1:K);

    % mIdx=6 is also good
    mIdx=5; as=1; k=1;
    mc=+inf;
    [mc,otherVar  ]=minAndName(nameList   ,mIdx,mc);
    [mc,pnpgd55Var]=minAndName(pnpgd55List,mIdx,mc);
    [mc,pnpg55Var ]=minAndName(pnpg55List ,mIdx,mc);
    [mc,pnpgdVar  ]=minAndName(pnpgdList  ,mIdx,mc);
    [mc,pnpgVar   ]=minAndName(pnpgList   ,mIdx,mc);
    [mc,tfocsVar  ]=minAndName(tfocsList  ,mIdx,mc);
    [mc,testVar   ]=minAndName(testList   ,mIdx,mc);
    [mc,spiralVar   ]=minAndName(spiralList   ,mIdx,mc);

    %compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgd55Var{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpg55Var{:});
    %compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgdVar{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),pnpgVar{:});
    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),testVar{:},spiralVar{:},tfocsVar{:});

    %compare({'cost'},@(y,varargin)semilogy((y-mc)/mc,varargin{:}),pnpgd55Var{:});
    compare({'cost'},@(y,varargin)semilogy((y-mc)/mc,varargin{:}),pnpg55Var{:});
    %compare({'cost'},@(y,varargin)semilogy((y-mc)/mc,varargin{:}),pnpgdVar{:});
    compare({'cost'},@(y,varargin)semilogy((y-mc)/mc,varargin{:}),pnpgVar{:});

%   compare({'innerItr'},@plot,varList{1:end/2});
%   compare({'innerItr'},@plot,varList{end/2:end});
%   return

    fields={'RMSE','time','cost'};
    forSave=addTrace(      pnpg_ {mIdx,as,k},     [],fields,mc); %  1- 4
    forSave=addTrace(     pnpg_d {mIdx,as,k},forSave,fields,mc); %  5- 8
    forSave=addTrace(      spiral{mIdx,as,k},forSave,fields,mc); %  9-12
    forSave=addTrace( pnpg_nInf  {mIdx,as,k},forSave,fields,mc); % 13-16
    forSave=addTrace(   pnpg_n0  {mIdx,as,k},forSave,fields,mc); % 17-20
    forSave=addTrace(tfocs_200_m6{mIdx,as,k},forSave,fields,mc); % 21-24
    forSave=addTrace(    pnpgA0  {mIdx,as,k},forSave,fields,mc); % 25-28
    forSave=addTrace(  pnpgG5A0  {mIdx,as,k},forSave,fields,mc); % 29-32
    forSave=addTrace(  pnpgG5Aq  {mIdx,as,k},forSave,fields,mc); % 33-26
    forSave=addTrace(  pnpgGfA0  {mIdx,as,k},forSave,fields,mc); % 37-40
    forSave=addTrace(  pnpgGfAq  {mIdx,as,k},forSave,fields,mc); % 41-44
    save('cost_itrPET.data','forSave','-ascii');

    fields_={'RMSE','time','cost'};
    forSave=addTrace(     pnpg_55{mIdx,as,k},     [],fields_,mc); %  1- 4
    forSave=addTrace(      cpdwt1{mIdx,as,k},forSave,fields_,mc); %  5- 8
    forSave=addTrace(      cpdwt2{mIdx,as,k},forSave,fields_,mc); %  9-12
    forSave=addTrace(tfocs_200_m9{mIdx,as,k},forSave,fields_,mc); % 13-16
    forSave=addTrace(       vmila{mIdx,as,k},forSave,fields_,mc); % 17-20
    save('cost_itrPET_1.data','forSave','-ascii');
    paperDir='~/research/myPaper/asilomar2014/';
    system(['mv cost_itrPET.data cost_itrPET_1.data ' paperDir]);

    return;

    mincost=reshape(forSave(:,[4,8,12,16,20,24]),[],1); 
    mincost=min(mincost(mincost~=0));

    figure; semilogy(forSave(:,5),'r'); hold on;
    semilogy(forSave(:,13),'b');
    semilogy(forSave(:,17),'k');
    %semilogy(forSave(:,9),'g');
    title('step size versus number of iterations');
    legend('pnpg','npg nInf','pnpg n0');

    figure;
    semilogy(forSave(:, 3),forSave(:, 4)-mincost,'r'); hold on;
    semilogy(forSave(:, 7),forSave(:, 8)-mincost,'g');
    semilogy(forSave(:,11),forSave(:,12)-mincost,'b');
    semilogy(forSave(:,15),forSave(:,16)-mincost,'k');
    semilogy(forSave(:,19),forSave(:,20)-mincost,'c');
    semilogy(forSave(:,23),forSave(:,24)-mincost,'k--');
    legend('npg n4','pnpg n4','spiral','pnpg nInf','pnpg n0','tfocsAT');
    hold on;

    keyboard

    idx=min(find(forSave(:,10)<1e-6));
    plot(forSave(idx,9),forSave(idx,7)-mincost,'bo');
    xxx=idx;
    idx=min(find(forSave(10:end,14)<1e-6))+10;
    plot(forSave(idx,13),forSave(idx,11)-mincost,'k*');
    xxx=[xxx;idx];  xxx=xxx(:)';
    save('cost_itrPETidx.data','xxx','-ascii');

    figure;
    semilogy(forSave(:, 3),forSave(:, 2),'r'); hold on;
    semilogy(forSave(:, 7),forSave(:, 6),'g');
    semilogy(forSave(:,11),forSave(:,10),'b');
    semilogy(forSave(:,15),forSave(:,14),'k');
    semilogy(forSave(:,19),forSave(:,18),'c');
    semilogy(forSave(:,23),forSave(:,22),'k--');
    legend('npg n4','pnpg n4','spiral','pnpg nInf','pnpg n0','tfocsAT');

    nn=128;
    xtrue = read_zubal_emis('nx', nn, 'ny', nn);
    % attenuation map
    mumap = read_zubal_attn('nx', nn, 'ny', nn);
    imwrite(xtrue/max(xtrue(:)),'pet.png');
    imwrite(mumap/max(mumap(:)),'mumap.png');

    idx=5;
    fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
    fprintf(' PNPGc: %g%%\n',pnpg_d{idx}.RMSE(end)*100);
    fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
    fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
    img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
    img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_pet.png')
    img=showImgMask(pnpg_d{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'PNPGc_pet.eps','psc2'); imwrite(img/max(xtrue(:)), 'PNPGc_pet.png')
    img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet.png')
    img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet.png')

%   idx=4;
%   fprintf('  PNPG: %g%%\n', pnpg_{idx}.RMSE(end)*100);
%   fprintf(' PNPGc: %g%%\n',pnpg_d{idx}.RMSE(end)*100);
%   fprintf('SPIRAL: %g%%\n',spiral{idx}.RMSE(end)*100);
%   fprintf('   FBP: (%g%%, %g%%)\n',   fbp{idx}.RMSE(end)*100,rmseTruncate(  fbp{idx},pnpg_{idx}.opt.trueX)*100);
%   img=pnpg_{idx}.x; mask=pnpg_{idx}.opt.mask;
%   img=showImgMask( pnpg_{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'PNPG_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),  'PNPG_pet2.png')
%   img=showImgMask(pnpg_d{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf, 'PNPGc_pet2.eps','psc2'); imwrite(img/max(xtrue(:)), 'PNPGc_pet2.png')
%   img=showImgMask(spiral{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'SPIRAL_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),'SPIRAL_pet2.png')
%   img=showImgMask(   fbp{idx}.x,mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,   'FBP_pet2.eps','psc2'); imwrite(img/max(xtrue(:)),   'FBP_pet2.png')

    paperDir='~/research/myPaper/asilomar2014/';
    decide=input(sprintf('start to copy to %s [y/N]?',paperDir),'s');
    if strcmpi(decide,'y')
        system(['mv varyCntPET.data cost_itrPET.data *_pet.png ' paperDir]);
    end
    system('rm *_pet.png *_pet.eps *_pet2.eps *_pet2.png');
    close all;

    figure;
    loglog(count,meanOverK(      fbp,'RMSE'),'b-o'); hold on;
    loglog(count,meanOverK(    pnpg_,'RMSE'),'r-*'); hold on;
    loglog(count,meanOverK(   spiral,'RMSE'),'k*-.');
    loglog(count,meanOverK(  pnpg_n0,'RMSE'),'c>-');
    loglog(count,meanOverK(pnpg_nInf,'RMSE'),'gs-');
    loglog(count,meanOverK(   pnpg_d,'RMSE'),'bp-.');
    loglog(count,meanOverK(    tfocs,'RMSE'),'kh-');
    legend('fbp','pnpg','spiral','pnpg\_n0','pnpg\_nInf','pnpg\_d','tfocs');

    figure;
    loglog(count,meanOverK(    pnpg_,'time'),'r-*'); hold on;
    loglog(count,meanOverK(   spiral,'time'),'k*-.');
    loglog(count,meanOverK(  pnpg_n0,'time'),'c>-');
    loglog(count,meanOverK(pnpg_nInf,'time'),'gs-');
    loglog(count,meanOverK(   pnpg_d,'time'),'bp-.');
    loglog(count,meanOverK(    tfocs,'time'),'kh-');
    legend('pnpg','spiral',' pnpg\_n0','pnpg\_nInf','pnpg\_d','tfocs');

    % time cost RMSE
    forSave=[count(:),meanOverK(   fbp,'RMSE'),...
        meanOverK(    pnpg_),...
        meanOverK(   spiral),...
        meanOverK(  pnpg_n0),...
        meanOverK(pnpg_nInf),...
        meanOverK(   pnpg_d),...
        meanOverK(    tfocs),...
        ];
    save('varyCntPET.data','forSave','-ascii');

    keyboard

case 'fullplot'
    filename = [mfilename '.mat'];
    load(filename);

    k=1;
    aa =(3:-0.5:-6);
    for i=1:length(count)
        pnpgContRMSE  {i,k} = [  pnpgFull{i,k}.contRMSE(:);  pnpgFull{i,k}.RMSE(end)]; out=pnpgContRMSE{i,k};
        fprintf('i=%d, good a = 1e%g PNPG\n',i,max((aa(out==min(out)))));
        fpnpgContRMSE {i,k} = [ fpnpgFull{i,k}.contRMSE(:); fpnpgFull{i,k}.RMSE(end)]; out=fpnpgContRMSE{i,k};
        fprintf('i=%d, good a = 1e%g FPNPG\n',i,max((aa(out==min(out)))));
        spiralContRMSE {i,k} = Cell.getField(spiralFull(i,:,k),'RMSE'); out=fpnpgContRMSE{i,k};
        fprintf('i=%d, good a = 1e%g SPIRAL\n',i,max((aa(out==min(out)))));
    end

    for i=1:length(count)
        figure;
        semilogy(aa(1:length(pnpgContRMSE{i})),pnpgContRMSE{i},'r-*'); hold on;
        semilogy(aa(1:length(fpnpgContRMSE{i})),fpnpgContRMSE{i},'g-o');
        semilogy(aa(1:length(spiralContRMSE{i})),spiralContRMSE{i},'b-s');
        title(num2str(i));
        legend('PNPG','FPNPG','SPIRAL');
        aaa(i)=min(pnpgContRMSE{i});
        bbb(i)=min(fpnpgContRMSE{i});
        ccc(i)=min(spiralContRMSE{i});
    end
    figure; semilogy(aaa,'r-*'); hold on;
    semilogy(bbb,'g-o');
    semilogy(ccc,'b-s');
    title('rmse vs count');
    legend('PNPG','FPNPG','SPIRAL');
end


ept=[];
pnpg_d_06=[];
pnpg_d55_06=[];
<<<<<<< HEAD
pnpg_adp=[];
pnpgG5A055=[];
pnpg_d_adp=[];
pnpg_1_3=[];
=======
ept=[];
pnpg_nInf88=[];
pnpgGfA088=[];
pnpgA0_d88=[];
pnpg_nInf_d88=[];
pnpgGfA0_d88 =[];
pnpg_d88=[];
pnpgG5A088=[];
pnpgGfAq88=[];
pnpg_n088=[];
pnpgG5Aq88=[];
pnpgGfAq_d88 =[];
pnpg_n0_d88=[];
pnpgG5Aq_d88=[];
spiral_m9=[];
pnpgA088=[];
pnpgG5A0_d88=[];
tfocs_200_m12=[];
tfocs_200_m6=[];
tfocs_200_m9=[];
>>>>>>> refs/remotes/origin/switch_to_object_proximal

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
    tt=getfield(method,fields{1});
    if(iscell(tt) && length(tt)==1)
        tt=tt{1};
    end
    itr=linspace(1,length(tt),len);
    if(~exist('fields','var'))
        fields={'time','cost','RMSE'};
    end
    for i=1:length(fields);
        tt=getfield(method,fields{i});
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
