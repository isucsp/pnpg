function PET_DWT_CP_Ex(op)
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
    %OPT.mask=[];
    OPT.outLevel=1;
    OPT.maxItr=2e3; OPT.thresh=1e-9; OPT.debugLevel=2; OPT.noiseType='poisson';
    OPT.maxItr=2e4; OPT.thresh=1e-9; OPT.debugLevel=2; OPT.noiseType='poisson';
    C.exact=true; C.val=@(x)0; C.prox=@(x,u)max(0,x);

    count = [1e4 1e5 1e6 1e7 1e8 1e9];
    K=1;

    as = [ 0.5, 0.5,0.5, 0.5, 0.5,   1];
    a  = [-0.5,   0,  0, 0.5, 0.5, 0.5];
    atv= [-0.5,-0.5,  0,   0, 0.5,   1];

    for k=1:K
    for i=5; %[1 2 3 4 6]
        j=1;
        [y,Phi,Phit,Psi,Psit,fbpfunc,OPT]=loadPET(count(i),OPT,k*100+i);
        NLL=@(x) Utils.poissonModel(x,Phi,Phit,y,OPT.bb);
        PROXOPT=[];
        PROXOPT.Lip=@(u)u^2; PROXOPT.initStep='fixed';
        PROXOPT.adaptiveStep=false; PROXOPT.backtracking=false;
        proximal=sparseProximal(Psi,Psit,C.prox,'pnpg',PROXOPT);
        proxOpt=PROXOPT;  proxOpt.dualGap=true;
        proximal_dualInnerCrit=sparseProximal(Psi,Psit,C.prox,'pnpg',proxOpt);

        fbp{i,1,k}.x=fbpfunc(y);
        fbp{i,1,k}.RMSE=sqrNorm(maskFunc(fbp{i,1,k}.x,OPT.mask)-OPT.trueX)/sqrNorm(OPT.trueX);

        fprintf('fbp RMSE=%f\n',fbp{i,1,k}.RMSE);
        fprintf('min=%d, max=%d, mean=%d\n',min(y(y>0)),max(y(y>0)),mean(y(y>0)));
        u_max=1;
        OPT.u = u_max*10.^a(i);

        initSig=C.prox(maskFunc(fbp{i,1,k}.x,OPT.mask));

        fprintf('%s, i=%d, j=%d, k=%d\n','PET Example',i,j,k);

        opt=OPT; opt.thresh=opt.thresh/100;   %opt.maxItr=3e3;
        L=1e5;                                %opt.xxx=pnpg_{i}.cost(end);
        sigma= [0,   0,   0,   0,1e-4];
        sigma1=[0,1e-3,1e-2,1e-1,1e-4];
        tau=   [1,   1,   1,   1,1e-4];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt25{i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);
        mysave
        return;

        opt=OPT; opt.thresh=opt.thresh/100;    
        L=1e5;                                 
        sigma= [0,   0,   0,   0,1e-6];
        sigma1=[0,1e-3,1e-2,1e-1,1e+1];
        tau=   [1,   1,   1,   1,1e-3];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt21{i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;    
        L=1e5;                                 
        sigma= [0,   0,   0,   0,1e-6];
        sigma1=[0,1e-3,1e-2,1e-1,1e-1];
        tau=   [1,   1,   1,   1,1e-3];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt22{i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;    
        L=1e5;                                 
        sigma= [0,   0,   0,   0,1e-6];
        sigma1=[0,1e-3,1e-2,1e-1,1e-3];
        tau=   [1,   1,   1,   1,1e-3];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt23{i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;    
        L=1e5;                                 
        sigma= [0,   0,   0,   0,1e-5];
        sigma1=[0,1e-3,1e-2,1e-1,1e-0];
        tau=   [1,   1,   1,   1,1e-3];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt24{i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.thresh=opt.thresh/100;  % opt.maxItr=3e3;
        L=1e5;                               % opt.xxx=pnpg_{i}.cost(end);
        sigma= [0,   0,   0,   0,10^-6];
        sigma1=[0,1e-3,1e-2,1e-1,    1];
        tau=   [1,   1,   1,   1,10^-3];
        opt.sigma=[sigma(i),sigma1(i)]; opt.tau=1/L/opt.sigma(1)*tau(i);
        cpdwt2 {i,j,k}=CP_DWT(Phi,Phit,y,3,Psi,Psit,C,initSig,opt);

        opt=OPT; opt.maxInnerItr=1e3;
        pnpg_ {i,j,k}=pnpg(NLL,proximal,initSig,opt);
        mysave

    end
    end

case lower('plot')
    filename = [mfilename '.mat'];
    load(filename);
    dup=load(filename);
    fprintf('PET Poisson l1 example\n');

    count = [1e4 1e5 1e6 1e7 1e8 1e9];

    cpdwt2x={
        'pnpg_',
        'cpdwt2',
        'cpdwt21',
        'cpdwt22',
        'cpdwt23',
        'cpdwt24',
        'cpdwt25',
        };

    K = 1;
    mIdx=5; as=1; k=1;
    mc=+inf;
    [mc,cpdwt2xVar]=minAndName(dup,cpdwt2x ,mIdx,mc);

    compare({'time','cost'},@(x,y,varargin)semilogy(x,(y-mc)/mc,varargin{:}),cpdwt2xVar{:});
    compare({'cost'},@(y,varargin)semilogy((y-mc)/mc,varargin{:}),cpdwt2xVar{:});

    fields={'RMSE','time','cost'};
    forSave=addTrace(pnpg_  {mIdx,as,k},     [],fields,mc); %  1- 4
    forSave=addTrace(cpdwt2 {mIdx,as,k},forSave,fields,mc); %  5- 8
    forSave=addTrace(cpdwt21{mIdx,as,k},forSave,fields,mc); %  9-12
    forSave=addTrace(cpdwt22{mIdx,as,k},forSave,fields,mc); % 13-16
    forSave=addTrace(cpdwt23{mIdx,as,k},forSave,fields,mc); % 17-20
    forSave=addTrace(cpdwt24{mIdx,as,k},forSave,fields,mc); % 21-24
    forSave=addTrace(cpdwt25{mIdx,as,k},forSave,fields,mc); % 25-28
    save('cpdwt2x.data','forSave','-ascii');
    paperDir='~/research/myPaper/asilomar2014/';
    system(['mv cpdwt2x.data ' paperDir]);
end

end

function [mc,varList] = minAndName(dup,nameList,i,mc)
    if(~exist('mc','var'))
        mc=+inf;
    end
    for ii=1:length(nameList)
        %a=eval(nameList{ii});
        a=dup.(strtrim(nameList{ii}));
        mc=min(mc,min(a{i}.cost(:)));
        a{i}.name=nameList{ii};
        varList{ii}=a{i};
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
    if(iscell(tt) && length(tt)==1)
        tt=tt{1};
    end
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
