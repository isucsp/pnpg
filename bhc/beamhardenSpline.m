function out = beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt)
%beamharden    beamharden effect correct method
%   out = beamharden(***)
%   Phi         The projection matrix implementation function handle
%   Phit        Transpose of Phi
%   Psi         Inverse wavelet transform matrix from wavelet coefficients
%               to image.
%   Psit        Transpose of Psi, need to have ΨΨ'=I
%   y           Log scale of Beamhardening measurement y=-log(I^{mea}/I_0)
%   xInit       Initial value for the algorithm
%   opt         Structure for the configuration of this algorithm (refer to
%               the code for detail)
%
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%
%   v_0.5:      use the Poisson measurement model as the default;
%               add GML model for the mass attenuation spectrum estimation
%   v_0.4:      use spline as the basis functions, make it more configurable
%   v_0.3:      add the option for reconstruction with known Ie
%   v_0.2:      add alphaDif to output;
%               add t[123] to output;
%
%   todo:       record the # of steps for the line search
%               make sure to add 1/2 to the likelihood
%               Try by have less number of sampling points.
%               use annihilating filter to do Ie estimation.
%               use cpu version of operators
%               optimize the form of Phi[t]Func51.m in subfuction
%

% for alpha step
% 'NPG'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
if(~isfield(opt,'noiseType')) opt.noiseType='Poisson'; end
if(~isfield(opt,'alphaStep')) opt.alphaStep='NPG'; end
if(~isfield(opt,'skipAlpha')) opt.skipAlpha=false; end
if(~isfield(opt,'maxAlphaSteps')) opt.maxAlphaSteps=1; end

if(~isfield(opt,'continuation')) opt.continuation=false; end
if(~isfield(opt,'contShrnk')) opt.contShrnk=0.98; end
if(~isfield(opt,'contCrtrn')) opt.contCrtrn=1e-4; end
if(~isfield(opt,'contEta')) opt.contEta=1e-2; end
if(~isfield(opt,'contGamma')) opt.contGamma=1e4; end

% default to not use nonnegative constraints.
if(~isfield(opt,'nu')) opt.nu=0; end
if(~isfield(opt,'u')) opt.u=1e-3; end
if(~isfield(opt,'restart')) opt.restart=true; end

if(~isfield(opt,'uMode')) opt.uMode='abs'; end
if(~isfield(opt,'a')) opt.a=1e-6; end

% for Ie step
if(~isfield(opt,'IeStep')) opt.IeStep='ActiveSet'; end
if(~isfield(opt,'skipIe')) opt.skipIe=false; end
if(~isfield(opt,'maxIeSteps')) opt.maxIeSteps=20; end
if(~isfield(opt,'etaDifCost')) opt.etaDifCost=1e-1; end
if(~isfield(opt,'spectBasis')) opt.spectBasis='dis'; end
if(~isfield(opt,'CenterB')) opt.CenterB=false; end

% The range for mass attenuation coeff is 1e-2 to 1e4 cm^2/g
if(~isfield(opt,'sampleMode')) opt.sampleMode='logspan'; end
if(~isfield(opt,'logspan')) opt.logspan=3; end
if(~isfield(opt,'K')) opt.K=2; end  % number of materials
if(~isfield(opt,'E')) opt.E=17; end

% common
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'initStep')) opt.initStep='BB'; end
if(~isfield(opt,'preSteps')) opt.preSteps=2; end
if(~isfield(opt,'errorType')) opt.errorType=0; end
% the higher, the more information. Set to 0 to turn off.
if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
if(~isfield(opt,'verbose')) opt.verbose=100; end
% Threshold for relative difference between two consecutive α
if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
if(~isfield(opt,'minItr')) opt.minItr=10; end   % currently not used
if(~isfield(opt,'saveAnimate')) opt.saveAnimate=false; end

% for NCG_PR
if(~isfield(opt,'muLustig')) opt.muLustig=1e-13; end

% use GML model to estimate Ie in after getting the estimation of α
if(~isfield(opt,'estIe')) opt.estIe=false; end

Imea=exp(-y); alpha=xInit(:); Ie=zeros(opt.E,1);

if(isfield(opt,'trueAlpha'))
    switch opt.errorType
        case 0
            trueAlpha = opt.trueAlpha/pNorm(opt.trueAlpha);
            computError= @(xxx) 1-(innerProd(xxx,trueAlpha)^2)/sqrNorm(xxx);
        case 1
            trueAlphaNorm=sqrNorm(opt.trueAlpha);
            if(trueAlphaNorm==0) trueAlphaNorm=eps; end
            computError = @(xxx) sqrNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
        case 2
            trueAlphaNorm=pNorm(opt.trueAlpha);
            if(trueAlphaNorm==0) trueAlphaNorm=eps; end
            computError = @(xxx) pNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
    end
end

if(opt.debugLevel>=3) figCost=1000; figure(figCost); end
if(opt.debugLevel>=4) figRes=1001; figure(figRes); end
if(opt.debugLevel>=5 && ~opt.skipIe) figIe=1003; figure(figIe); end
if(opt.debugLevel>=6) figAlpha=1002; figure(figAlpha); end

switch lower(opt.sampleMode)
    case 'assigned'
        Ie=zeros(length(opt.kappa),1);
        temp=opt.kappa;
        temp1=abs(temp-1);
        temp2=find(temp1==min(temp1));
        Ie(temp2-1:temp2+1)=1/3;
    case 'logspan'
        temp=logspace(-floor(opt.E/2)/(opt.E-1)*opt.logspan,...
            floor(opt.E/2-0.5)/(opt.E-1)*opt.logspan,opt.E);
        opt.E=length(temp);
        Ie=zeros(opt.E,1); Ie(floor(opt.E/2)+1)=1;
end

kappa=temp(:); q=kappa(2)/kappa(1);
temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
polymodel = Spline(opt.spectBasis,temp);
polymodel.setPlot(opt.kappa,opt.iota,opt.epsilon);
polyIout = polymodel.polyIout;

Ie=Ie/polyIout(0,Ie);

[opt.upkappa,opt.upiota]=getUpiota(opt.epsilon,opt.kappa,opt.iota);

% find the best intial Ie ends
if(isfield(opt,'Ie')) Ie=opt.Ie(:);
else
    if(opt.skipIe)  % it is better to use dis or b-1 spline
        Ie=interp1(log(opt.upkappa), opt.upiota ,log(kappa(:)),'spline');
        % there will be some points interplated negative and need to be removed
        Ie=max(Ie,0);
    end
end

switch lower(opt.alphaStep)
    case lower('NCG_PR')
        alphaStep = NCG_PR(3,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        if(isfield(opt,'muHuber'))
            fprintf('use huber approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) Utils.huber(aaa,opt.muHuber,Psi,Psit);
        end
        if(isfield(opt,'muLustig'))
            fprintf('use lustig approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) Utils.lustigL1(aaa,opt.muLustig,Psi,Psit);
        end
    case lower('NPGs')
        alphaStep=NPGs(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
    case lower('NPG')
        alphaStep=NPG (1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
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

switch lower(opt.IeStep)
    case lower('ActiveSet')
        IeStep = ActiveSet(B,b,Ie,opt.maxIeSteps,opt.stepShrnk);
    case lower('NPG')
        IeStep = NPG_ind(Ie,true,[],[],opt.maxIeSteps,opt.stepShrnk,opt.thresh);
end

switch lower(opt.noiseType)
    case lower('Gaussian')
        alphaStepFunc = @(III,aaa,pIout) gaussLAlpha(Imea,III,aaa,Phi,Phit,pIout,IeStep);
        IeStepFunc = @(A,III) gaussLI(Imea,A,III);
    case lower('Poisson')
        alphaStepFunc = @(III,aaa,pIout) poissLAlpha(Imea,III,aaa,Phi,Phit,pIout);
        IeStepFunc = @(A,III) poissLI(Imea,A,III);
end
alphaStep.fArray{1} = @(aaa) alphaStepFunc(Ie,aaa,polyIout);

if(strcmpi(opt.uMode,'relative'))
    PsitPhitz=Psit(Phit(y));
    PsitPhit1=Psit(Phit(ones(length(y),1)));
    temp=[]; [temp(1),temp(2)]=polyIout(0,Ie);
    opt.u=opt.a*max(abs(PsitPhitz+PsitPhit1*log(temp(1))))*temp(2)/temp(1);
end
if(opt.continuation)
    PsitPhitz=Psit(Phit(y));
    PsitPhit1=Psit(Phit(ones(length(y),1)));
    temp=[]; [temp(1),temp(2)]=polyIout(0,Ie);
    alphaStep.u=0.1*max(abs(PsitPhitz+PsitPhit1*log(temp(1))))*temp(2)/temp(1);
    % This makes sure that the regulator settles after 300 iterations
    alphaStep.u = min(alphaStep.u,opt.u*1000);

    keyboard
    contIdx=1;
    [~,g]=alphaStep.fArray{1}(alpha);
    inf_psit_grad=pNorm(Psit(g),inf);
    alphaStep.u = opt.contEta*inf_psit_grad;
    alphaStep.u = min(alphaStep.u,opt.u*opt.contGamma);
    alphaStep.u = max(alphaStep.u,opt.u);
    if(alphaStep.u*opt.contShrnk<=opt.u)
        opt.continuation=false;
        alphaStep.u=opt.u;
    end
    clear('g');
else alphaStep.u=opt.u;
    fprintf('opt.u=%g\n',opt.u);
end

disp(['use initial sparsity regulator u:' num2str(alphaStep.u)]);

if(any(strcmp(properties(alphaStep),'restart')))
    if(~opt.restart) alphaStep.restart=-1; end
else
    opt.restart=false;
end

if(any(strcmp(properties(alphaStep),'adaptiveStep'))...
        && isfield(opt,'adaptiveStep'))
    % use the default of alphaStep
    alphaStep.adaptiveStep=opt.adaptiveStep;
end
if(any(strcmp(properties(alphaStep),'cumuTol'))...
        && isfield(opt,'cumuTol'))
    alphaStep.cumuTol=opt.cumuTol;
end

if(strcmpi(opt.initStep,'fixed'))
    alphaStep.stepSizeInit(opt.initStep,opt.L);
else alphaStep.stepSizeInit(opt.initStep);
end

if(any(strcmp(properties(alphaStep),'innerSearch')))
    collectInnerSearch=true;
else
    collectInnerSearch=false;
end

if(any(strcmp(properties(alphaStep),'debug')))
    collectDebug=true; out.debug={};
else
    collectDebug=false;
end

if(any(strcmp(properties(alphaStep),'preSteps')))
    alphaStep.preSteps=opt.preSteps;
end

tic; p=0; strlen=0; convThresh=0;
while( ~(opt.skipAlpha && opt.skipIe) )
    if(opt.saveAnimate && (mod(p,10)==0 || p<10))
        img=showImgMask(alpha,opt.mask);
        imwrite(img/max(img(:)),sprintf('animate_%03d.png',p),'png');
        [~,~,kkkkk,temp]=polymodel.plotSpectrum(Ie);
        if(p==0) IeAnimate=kkkkk(:); end
        IeAnimate=[IeAnimate, temp(:)];
        save('IeAnimate.data','IeAnimate','-ascii');
    end

    p=p+1; str=sprintf('p=%-4d',p);
    
    % start optimize over alpha
    if(~opt.skipAlpha) 
        alphaStep.fArray{1} = @(aaa) alphaStepFunc(Ie,aaa,polyIout);
        alphaStep.main();
        
        out.fVal(p,:) = (alphaStep.fVal(:))';
        out.cost(p) = alphaStep.cost;
        if(~opt.skipIe) IeStep.cost=alphaStep.fVal(1); end

        out.alphaSearch(p) = alphaStep.ppp;
        out.stepSize(p) = alphaStep.stepSize;
        if(alphaStep.restart>=0) out.restart(p)=alphaStep.restart; end
        if(collectInnerSearch) out.innerSearch(p)=alphaStep.innerSearch; end;
        if(collectDebug && ~isempty(alphaStep.debug))
            out.debug{size(out.debug,1)+1,1}=p;
            out.debug{size(out.debug,1),2}=alphaStep.debug;
        end;
        if(opt.debugLevel>1)
            out.BB(p,1)=alphaStep.stepSizeInit('BB');
            out.BB(p,2)=alphaStep.stepSizeInit('hessian');
            % alphaStep.stepSizeInit('hessian',alpha);
        end

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
        if(strcmpi(opt.noiseType,'Gaussian'))
            str=sprintf([str ' zmf=(%g,%g)'], IeStep.zmf(1), IeStep.zmf(2));
            out.zmf(p,:)=IeStep.zmf(:)';
        end
        if(p>1)
            str=sprintf([str ' difCost=%g'], out.difCost(p));
        end
    end
    % end optimizing over alpha
    
    %if(out.delta<=1e-4) maxPP=5; end
    if(~opt.skipIe)
        % update the object fuction w.r.t. Ie
        if(strcmpi(opt.noiseType,'Gaussian') && strcmpi(opt.IeStep,'NPG'))
            if(max(IeStep.zmf(:))>=1 || p<20)  IeStep.maxItr=opt.maxIeSteps*10;
            else IeStep.maxItr=opt.maxIeSteps; end
        end

        A = polyIout(Phi(alpha),[]);
        IeStep.func = @(III) IeStepFunc(A,III);
        if(p>1) IeStep.thresh=opt.etaDifCost*out.difCost(p); end
        [~,IeStepCnt]=IeStep.main();

        alphaStep.cost=IeStep.cost+opt.u*alphaStep.fVal(3);
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
            case lower('NPG')
                out.IeSearch(p) = IeStep.ppp;
                str=sprintf([str ' IeSearch=%d'],out.IeSearch(p));
                str=sprintf([str ' #IeStep=%d'],IeStepCnt);
        end

        out.difIe(p)=norm(IeStep.Ie(:)-Ie(:))/norm(Ie);
        
        if(p>1) out.difllI(p)=abs(out.llI(p)-out.llI(p-1))/out.llI(p); end
        Ie = IeStep.Ie;
        str=sprintf([str ' difIe=%g'], out.difIe(p));
        if(p>1) str=sprintf([str ' difllI=%g'],out.difllI(p)); end
    end

    if(opt.debugLevel>1)
        if(p>1 && opt.debugLevel>=3)
            set(0,'CurrentFigure',figCost);
            if(isfield(opt,'trueAlpha')) subplot(2,1,1); end
            if(out.cost(p)>0)
                semilogy(p-1:p,out.cost(p-1:p),'k'); hold on;
                title(sprintf('cost(%d)=%g',p,out.cost(p)));
            end

            if(isfield(opt,'trueAlpha'))
                subplot(2,1,2);
                semilogy(p-1:p,out.RMSE(p-1:p)); hold on;
                title(sprintf('RMSE(%d)=%g',p,out.RMSE(p)));
            end
            drawnow;
        end
        if(~opt.skipIe && opt.debugLevel>=5)
            set(0,'CurrentFigure',figIe);
            polymodel.plotSpectrum(Ie);
            title(sprintf('int upiota d kappa = %g',polyIout(0,Ie)));
            drawnow;
        end
        if(length(out.fVal(p,:))>=1 && p>1 && opt.debugLevel>=4)
            set(0,'CurrentFigure',figRes);
            style={'r','g','b'};
            for i=1:length(out.fVal(p,:))
                subplot(3,1,i);
                semilogy(p-1:p,out.fVal(p-1:p,i),style{i}); hold on;
            end
            drawnow;
        end
        if(opt.debugLevel>=6)
            set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
            if(~isempty(IeStep.zmf))
                title(sprintf('zmf=(%g,%g)', IeStep.zmf(1), IeStep.zmf(2)))
            end
            drawnow;
        end
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
    end
    out.time(p)=toc;
    if( p>1 && out.difAlpha(p)<=opt.thresh && (alphaStep.u==opt.u(end)) )
        convThresh=convThresh+1;
    end
    if( p >= opt.maxItr || convThresh>2 )
        if(opt.debugLevel==0) fprintf('%s',str); end
        break;
    end
end

out.Ie=Ie; out.kappa=kappa;
out.alpha=alpha; out.p=p; out.opt = opt;
out.grad=alphaStep.grad;
out.date=datestr(now);
fprintf('\nTime used: %d, cost=%g',out.time(end),out.cost(end));
if(isfield(opt,'trueAlpha'))
        fprintf(', RMSE=%g\n',out.RMSE(end));
else
    fprintf('\n');
end

if(opt.estIe)
end

end

