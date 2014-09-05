function out = beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt)
%beamharden    beamharden effect correct method
%   out = beamharden(***)
%   Phi         The projection matrix implementation function handle
%   Phit        Transpose of Phi
%   Psi         Inverse wavelet transform matrix from wavelet coefficients
%               to image.
%   Psit        Transpose of Psi
%   y           Log scale of Beamhardening measurement y=-log(I^{mea}/I_0)
%   xInit       Initial value for the algorithm
%   opt         Structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%
%   v_0.4:      use spline as the basis functions, make it more configurable
%   v_0.3:      add the option for reconstruction with known Ie
%   v_0.2:      add llAlphaDif to output;
%               add t[123] to output;
%
%   todo:       record the # of steps for the line search
%               make sure to add 1/2 to the likelihood
%               Try by have less number of sampling points.
%               use annihilating filter to do Ie estimation.
%               use cpu version of operators
%               optimize the form of Phi[t]Func51.m in subfuction
%

%opt.a=-6.5;  % aArray=-6.8:0.2:-6.2;
if(~isfield(opt,'K')) opt.K=2; end
if(~isfield(opt,'E')) opt.E=17; end
% default to not use nonnegative constraints.
if(~isfield(opt,'nu')) opt.nu=0; end
if(~isfield(opt,'u')) opt.u=1e-3; end
if(~isfield(opt,'uMode')) opt.uMode='abs'; end
if(~isfield(opt,'a')) opt.a=1e-6; end
% the number of non-zeros in the transform domain
if(~isfield(opt,'rInit')) opt.rInit=5000; end
if(~isfield(opt,'spectBasis')) opt.spectBasis='dis'; end
if(~isfield(opt,'muLustig')) opt.muLustig=1e-13; end
if(~isfield(opt,'numCall')) opt.numCall=1; end
if(~isfield(opt,'logspan')) opt.logspan=3; end
if(~isfield(opt,'skipAlpha')) opt.skipAlpha=false; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'initStep')) opt.initStep='BB'; end
if(~isfield(opt,'skipIe')) opt.skipIe=false; end
% The range for mass attenuation coeff is 1e-2 to 1e4 cm^2/g
if(~isfield(opt,'muRange')) opt.muRange=[1e-2; 1e4]; end
if(~isfield(opt,'sampleMode')) opt.sampleMode='logspan'; end
if(~isfield(opt,'maxAlphaSteps')) opt.maxAlphaSteps=1; end
if(~isfield(opt,'maxIeSteps')) opt.maxIeSteps=20; end
if(~isfield(opt,'showImg')) opt.showImg=false; end
% the higher, the more information. Set to 0 to turn off.
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% 'FISTA_ADMM_NNL1'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
if(~isfield(opt,'alphaStep')) opt.alphaStep='FISTA_ADMM_NNL1'; end
if(~isfield(opt,'IeStep')) opt.IeStep='ActiveSet'; end
if(~isfield(opt,'continuation')) opt.continuation=false; end
if(~isfield(opt,'contShrnk')) opt.contShrnk=0.98; end
if(~isfield(opt,'contCrtrn')) opt.contCrtrn=1e-4; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-12; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'errorType')) opt.errorType=0; end
if(~isfield(opt,'saveAnimate')) opt.saveAnimate=false; end

Imea=exp(-y); alpha=xInit(:); Ie=zeros(opt.E,1);

if(isfield(opt,'trueAlpha'))
    trueAlphaNorm=norm(opt.trueAlpha);
    switch opt.errorType
        case 0
            trueAlpha = opt.trueAlpha/trueAlphaNorm;
            computError= @(xxx) 1-(xxx(:)'*trueAlpha/norm(xxx(:)))^2;
        case 1
            computError = @(xxx) (norm(xxx(:)-opt.trueAlpha)/trueAlphaNorm)^2;
    end
end

if(opt.debugLevel>=2) figCost=1000; figure(figCost); end
if(opt.debugLevel>=3) figRes=1001; figure(figRes); end
if(opt.debugLevel>=3 && ~opt.skipIe) figIe=1003; figure(figIe); end
if(opt.showImg && opt.debugLevel>=6) figAlpha=1002; figure(figAlpha); end

switch lower(opt.sampleMode)
    case 'uniform'
        temp=linspace(opt.muRange(1),opt.muRange(2),opt.E);
        Ie(floor(opt.E/2)-1:floor(opt.E/2)+1)=1/3;
    case 'exponential'
        temp=logspace(log10(opt.muRange(1)),log10(opt.muRange(2)),opt.E);
        temp1=abs(temp-1);
        Ie(temp1==min(temp1))=1;
    case 'assigned'
        Ie=zeros(length(opt.kappa),1);
        temp=opt.kappa;
        temp1=abs(temp-1);
        temp2=find(temp1==min(temp1));
        Ie(temp2-1:temp2+1)=1/3;
    case 'logspan'
        temp=logspace(-floor(opt.E/2)/(opt.E-1)*opt.logspan,...
            floor(opt.E/2-0.5)/(opt.E-1)*opt.logspan,opt.E);
        % q=(1+sqrt(5))/2; temp1=ceil(opt.logspan/2/log10(q));
        % temp=q.^(-temp1:temp1);
        % opt.E=2*temp1+1;
        opt.E=length(temp);
        Ie=zeros(opt.E,1); Ie(floor(opt.E/2)+1)=1;
end

kappa=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(opt.K-1)*eps);
temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
polymodel = Spline(opt.spectBasis,temp);
polymodel.setPlot(opt.trueKappa,opt.trueIota,opt.epsilon);
polyIout = polymodel.polyIout;

Ie=Ie/polyIout(0,Ie);

temp1 = [opt.epsilon(1);(opt.epsilon(1:end-1)+opt.epsilon(2:end))/2;opt.epsilon(end)];
temp2 = [opt.trueKappa(1);(opt.trueKappa(1:end-1)+opt.trueKappa(2:end))/2;opt.trueKappa(end)];
temp1 = temp1(2:end)-temp1(1:end-1);
temp2 = temp2(2:end)-temp2(1:end-1);
opt.trueUpiota=abs(opt.trueIota.*temp1./temp2);
opt.trueIota=opt.trueIota/(opt.trueIota'*temp1);
clear 'temp1' 'temp2'

% find the best intial Ie ends
if(isfield(opt,'Ie')) Ie=opt.Ie(:);
else
    if(opt.skipIe)  % it is better to use dis or b-1 spline
        Ie=interp1(opt.trueKappa, opt.trueUpiota,kappa(:),'spline');
        % there will be some points interplated negative and need to be removed
        Ie(Ie<0)=0;
    end
end

% find the best intial Ie starts
% R = polyIout(Phi(alpha),[]);
% for i=1:size(R,2)
%     temp(i) = var(y+log(R(:,i)),1);
% end
% idx = find(temp==min(temp));
% Ie = Ie*0;
% Ie(idx) = exp(-mean(y+log(R(:,idx))));
%max(Imea./(exp(-atten(Phi,alpha)*kappa')*Ie))

switch lower(opt.alphaStep)
    case lower('NCG_PR')
        alphaStep = NCG_PR(3,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        if(isfield(opt,'muLustig'))
            fprintf('use lustig approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) Utils.lustigL1(aaa,opt.muLustig,Psi,Psit);
        end
        if(isfield(opt,'muHuber'))
            fprintf('use huber approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) Utils.huber(aaa,opt.muHuber,Psi,Psit);
        end
    case {lower('SpaRSA')}
        alphaStep = Methods(2,alpha);
        alphaStep.coef(1:2) = [1; 1];
        alphaStep.fArray{2} = nonneg;
        alphaStep.maxStepNum = opt.maxAlphaSteps;
        alphaStep.stepShrnk = opt.stepShrnk;
        alphaStep.Psi = Psi;
        alphaStep.Psit = Psit;
        alphaStep.M = 5;
    case lower('FISTA_L1')
        alphaStep=FISTA_L1(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
    case lower('FISTA_ADMM_NNL1')
        alphaStep=FISTA_ADMM_NNL1(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
    case lower('IST_ADMM_NNL1')
        alphaStep=IST_ADMM_NNL1(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_NNL1')}
        alphaStep=ADMM_NNL1(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
    case {lower('ADMM_L1')}
        alphaStep = ADMM_L1(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
end
alphaStep.fArray{2} = @Utils.nonnegPen;
alphaStep.coef(1:2) = [1; opt.nu;];

B=eye(opt.E); b=zeros(opt.E,1);
if(isfield(opt,'CenterB') && opt.CenterB)
    if(~isfield(opt,'correctCenterB')) opt.correctCenterB=true; end
    temp=-eye(opt.E); temp(floor(opt.E/2)+1,:)=[]; temp(:,floor(opt.E/2)+1)=1;
    temp=temp/sqrt(2);
    B=[B; temp]; b=[b; zeros(opt.E-1,1)];
end
temp = polyIout(0,[]);
B=[B; -temp(:)'/norm(temp)]; b=[b; -1/norm(temp)];

switch lower(opt.IeStep)
    case lower('ActiveSet')
        IeStep = ActiveSet(B,b,Ie,opt.maxIeSteps,opt.stepShrnk);
    case lower('FISTA')
        IeStep = FISTA_Simplex(B,b,Ie,opt.maxIeSteps,opt.stepShrnk);
end

PsitPhitz=Psit(Phit(y));
PsitPhit1=Psit(Phit(ones(length(y),1)));
if(strcmpi(opt.uMode,'relative'))
    temp=[]; [temp(1),temp(2)]=polyIout(0,Ie);
    opt.u=opt.a*max(abs(PsitPhitz+PsitPhit1*log(temp(1))))*temp(2)/temp(1);
end
if(opt.continuation)
    temp=[]; [temp(1),temp(2)]=polyIout(0,Ie);
    alphaStep.u=0.1*max(abs(PsitPhitz+PsitPhit1*log(temp(1))))*temp(2)/temp(1);
    % This makes sure that the regulator settles after 300 iterations
    alphaStep.u = min(alphaStep.u,opt.u*1000);
else alphaStep.u=opt.u;
end
disp(['use initial sparsity regulator u:' num2str(alphaStep.u)]);

if(any(strcmp(properties(alphaStep),'restart')) ...
        && isfield(opt,'restart') && (~opt.restart))
    alphaStep.restart=-1;
end
if(any(strcmp(properties(alphaStep),'adaptiveStep'))...
        && isfield(opt,'adaptiveStep'))
    alphaStep.adaptiveStep=opt.adaptiveStep;
end
if(any(strcmp(properties(alphaStep),'cumuTol'))...
        && isfield(opt,'cumuTol'))
    alphaStep.cumuTol=opt.cumuTol;
end

alphaStep.fArray{1} = @(aaa) gaussLAlpha(Imea,Ie,aaa,Phi,Phit,polyIout,IeStep);
if(strcmpi(opt.initStep,'fixed'))
    alphaStep.stepSizeInit(opt.initStep,opt.L);
else alphaStep.stepSizeInit(opt.initStep);
end

tic; p=0; str=''; strlen=0; convThresh=0;
while( ~(opt.skipAlpha && opt.skipIe) )
    if(opt.saveAnimate && (mod(p,10)==0 || p<10))
        img=showImgMask(alpha,opt.mask);
        imwrite(img/max(img(:)),sprintf('animate_%03d.png',p),'png');
        [~,~,kkkkk,temp]=polymodel.plotSpectrum(Ie);
        if(p==0) IeAnimate=kkkkk(:); end
        IeAnimate=[IeAnimate, temp(:)];
        save('IeAnimate.data','IeAnimate','-ascii');
    end

    p=p+1;
    str=sprintf('p=%-4d',p);
    
    % start optimize over alpha
    if(~opt.skipAlpha)
        alphaStep.fArray{1} = @(aaa) gaussLAlpha(Imea,Ie,aaa,Phi,Phit,polyIout,IeStep);
        alphaStep.main();
        
        out.fVal(p,:) = (alphaStep.fVal(:))';
        out.cost(p) = alphaStep.cost;
        out.alphaSearch(p) = alphaStep.ppp;
        out.stepSize(p) = alphaStep.stepSize;
        if(alphaStep.restart>=0) out.restart(p)=alphaStep.restart; end

        out.difAlpha(p)=relativeDif(alphaStep.alpha,alpha);
        if(p>1) out.difCost(p)=abs(out.cost(p)-out.cost(p-1))/out.cost(p); end
        alpha = alphaStep.alpha;

        if(opt.continuation || strcmpi(opt.uMode,'relative'))
            out.uRecord(p,:)=[opt.u,alphaStep.u];
            str=sprintf([str ' u=%-6g'],alphaStep.u);
        end
        if(strcmpi(opt.uMode,'relative') && ~opt.skipIe)
            temp=[]; [temp(1),temp(2)]=polyIout(0,Ie);
            opt.u=opt.a*max(abs(PsitPhitz+PsitPhit1*log(temp(1))))*temp(2)/temp(1);
        end
        if(opt.continuation)
            %if(p>1 && opt.difCost(p)<opt.contCrtrn)
                alphaStep.u = max(alphaStep.u*opt.contShrnk,opt.u);
            %end
        else alphaStep.u=opt.u;
        end
        if(isfield(opt,'trueAlpha'))
            out.RMSE(p)=computError(alpha);
        end

        str=sprintf([str ' cost=%-6g RSE=%g',...
            ' difAlpha=%g aSearch=%d'],...
            out.cost(p),out.RMSE(p), out.difAlpha(p), ...
            alphaStep.ppp);
        str=sprintf([str ' zmf=(%g,%g)'], IeStep.zmf(1), IeStep.zmf(2));
        if(p>1)
            str=sprintf([str ' difCost=%g'], out.difCost(p));
        end
    end
    % end optimizing over alpha
    
    %if(out.delta<=1e-4) maxPP=5; end
    if(~opt.skipIe && ((~opt.skipAlpha && max(IeStep.zmf(:))<1) || (opt.skipAlpha)))
        % update the object fuction w.r.t. Ie
        A = polyIout(Phi(alpha),[]);
        IeStep.func = @(III) gaussLI(Imea,A,III);
        IeStep.main();

        if(isfield(opt,'CenterB') && opt.CenterB && opt.correctCenterB)
            temp=find(IeStep.Q((1:length(Ie)-1)+length(Ie)));
            if(~isempty(temp) && isContinuous(temp) && (any(temp==floor(length(Ie)/2)) || any(temp==floor(length(Ie)/2)+1)))
                temp=(mean(temp)-floor(length(Ie)/2)-0.5);
                q=kappa(2)/kappa(1);
                if temp>1
                    alphaStep.alpha=alphaStep.alpha*q;
                    IeStep.Ie=[IeStep.Ie(2:end); 0]*q;
                elseif  temp<-1
                    alphaStep.alpha=alphaStep.alpha/q;
                    IeStep.Ie=[0; IeStep.Ie(1:end-1)]/q;
                end
                if(abs(temp)>1)
                    alphaStep.reset();
                    IeStep.init();
                end
            end
        end

        out.llI(p) = IeStep.cost;
        switch lower(opt.IeStep)
            case lower('ActiveSet')
                out.IeSteps(p)= IeStep.stepNum;
                out.course{p} = IeStep.course;
                str=sprintf([str ' #steps=%d'],out.IeSteps(p));
                if(IeStep.converged)
                    str = [str ' cnvrgd'];
                else
                    str = [str ' uncnvrgd'];
                end
            case lower('FISTA')
                out.IeSearch(p) = IeStep.ppp;
                str=sprintf([str ' IeSearch=%d'],out.IeSearch(p));
        end
        out.difIe(p)=norm(IeStep.Ie(:)-Ie(:))/norm(Ie);
        if(p>1) out.difllI(p)=abs(out.llI(p)-out.llI(p-1))/out.llI(p); end
        Ie = IeStep.Ie;
        str=sprintf([str ' difIe=%g'], out.difIe(p));
        if(p>1) str=sprintf([str ' difllI=%g'],out.difllI(p)); end
    end

    if(p>1 && opt.debugLevel>=2)
        set(0,'CurrentFigure',figCost);
        if(isfield(opt,'trueAlpha')) subplot(2,1,1); end
        semilogy(p-1:p,out.cost(p-1:p),'k'); hold on;
        title(sprintf('cost(%d)=%g',p,out.cost(p)));

        if(isfield(opt,'trueAlpha'))
            subplot(2,1,2);
            semilogy(p-1:p,out.RMSE(p-1:p)); hold on;
            title(sprintf('RMSE(%d)=%g',p,out.RMSE(p)));
        end
        drawnow;
    end
    
    if(~opt.skipIe && opt.debugLevel>=3)
        set(0,'CurrentFigure',figIe);
        polymodel.plotSpectrum(Ie);
        title(sprintf('int upiota d kappa = %g',polyIout(0,Ie)));
        drawnow;
    end
    if(length(out.fVal(p,:))>=1 && p>1 && opt.debugLevel>=3)
        set(0,'CurrentFigure',figRes);
        style={'r','g','b'};
        for i=1:length(out.fVal(p,:))
            subplot(3,1,i);
            semilogy(p-1:p,out.fVal(p-1:p,i),style{i}); hold on;
        end
        drawnow;
    end
    if(opt.showImg && opt.debugLevel>=6)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        if(~isempty(IeStep.zmf))
            title(sprintf('zmf=(%g,%g)', IeStep.zmf(1), IeStep.zmf(2)))
        end
        drawnow;
    end
    
    if(opt.debugLevel>=1)
        if(alphaStep.warned || IeStep.warned)
            fprintf('%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
        if(mod(p,opt.verbose)==0)
            strlen=0; fprintf('\n');
        else strlen = length(str);
        end
        str='';
    end
    out.time(p)=toc;
    if(p>1 && out.difAlpha(p)<opt.thresh && (opt.skipIe || (out.difIe(p)<opt.thresh)) )
        convThresh=convThresh+1;
    end
    if(p >= opt.maxItr || convThresh>10)
        if(opt.debugLevel==0) fprintf('%s',str); end
        break;
    end
end

out.Ie=Ie; out.kappa=kappa; out.alpha=alpha; out.cpuTime=toc; out.p=p;
out.opt = opt;
out.Phi=Phi; out.Phit=Phit; out.Psi=Psi; out.Psit=Psit; out.y=y; out.xInit=xInit;
out.grad=alphaStep.grad;

fprintf('\n');

end

function out=isContinuous(x)
    y=x(2:end)-x(1:end-1);
    if(any(y~=1)) out=false; else out=true; end
end

