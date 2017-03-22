function PET_TV_CP_Ex(op)
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
    %for i=4:5
    %for i=[1,2,3,5,6]
    for i = 4;
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

        opt=OPT; opt.thresh=opt.thresh/100;       %opt.maxItr=3e3; opt.xxx=pnpg_{i}.cost(end);
        sigma =[ 1e2,1,1e-2, 10^-3,1e-6];
        sigma1=[   1,1,   1, 10^-3,1];
        tau=   [1e-1,1,   1, 10^-3,10^-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv25 {i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);
        mysave
        return

        opt=OPT;
        pnpg_ {i,j,k}=pnpg(NLL,proximal,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;
        sigma =[1e2,1,1e-2,1e-4,1e-6];
        sigma1=[1,1,1,1,1];
        tau=[1e-1,1,1,10^-2,10^-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv2 {i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;
        sigma =[ 1e2,1,1e-2, 10^-4,1e-6];
        sigma1=[   1,1,   1, 10^-1,1];
        tau=   [1e-1,1,   1, 10^-2,10^-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv21{i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;
        sigma =[ 1e2,1,1e-2, 10^-4,1e-6];
        sigma1=[   1,1,   1, 10^-2,1];
        tau=   [1e-1,1,   1, 10^-2,10^-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv22{i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;
        sigma =[ 1e2,1,1e-2, 10^-4,1e-6];
        sigma1=[   1,1,   1, 10^-3,1];
        tau=   [1e-1,1,   1, 10^-3,10^-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv23{i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;
        sigma =[ 1e2,1,1e-2, 10^-4,1e-6];
        sigma1=[   1,1,   1, 10^-4,1];
        tau=   [1e-1,1,   1, 10^-4,10^-3];
        opt.sigma=[sigma(i),sigma1(i),sigma1(i)]; opt.tau=tau(i);
        cptv24{i,j,k}=CP_TV(Phi,Phit,y,2,tvType,C,initSig,opt);

        mysave
    end
    end

case lower('plot')
    filename = [mfilename '.mat']; load(filename);
    fprintf('PET Poisson TV example for TSP\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    cptv2x={
        'pnpg_',
        'cptv2',
        'cptv21',
        'cptv22',
        'cptv23',
        'cptv24'};

    i=4;
    mc=+inf;
    [mc,cptv2xVar]=minAndName(cptv2x,i,mc);
    compare({'cost'},@(y,varargin)semilogy((y-mc)/mc,varargin{:}),cptv2xVar{:});

    fields_={'RMSE','time','cost'};
    forSave=addTrace(    pnpg_{i},     [],fields_,mc); %  1- 4
    forSave=addTrace(   cptv2 {i},forSave,fields_,mc); %  5- 8
    forSave=addTrace(   cptv21{i},forSave,fields_,mc); %  9-12
    forSave=addTrace(   cptv22{i},forSave,fields_,mc); % 13-16
    forSave=addTrace(   cptv23{i},forSave,fields_,mc); % 17-20
    forSave=addTrace(   cptv24{i},forSave,fields_,mc); % 21-24
    save('cptv2x.data','forSave','-ascii');
    paperDir='~/research/myPaper/asilomar2014/';
    system(['mv cptv2x.data ' paperDir]);
end

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

