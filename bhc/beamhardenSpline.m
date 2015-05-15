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
%               switch the order of alpha and I step;
%               use relative cost difference as convergence for I step
%               add GML model for the mass attenuation spectrum estimation
%   v_0.4:      use spline as the basis functions, make it more configurable
%   v_0.3:      add the option for reconstruction with known Ie
%   v_0.2:      add alphaDif to output;
%               add t[123] to output;
%

% for alpha step
% 'NPG'; %'SpaRSA'; %'NCG_PR'; %'ADMM_L1'; %
if(~isfield(opt,'noiseType')) opt.noiseType='Poisson'; end
if(~isfield(opt,'proximal' )) opt.proximal='wvltADMM'; end
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
if(~isfield(opt,'IeStep')) opt.IeStep='LBFGSB'; end
if(~isfield(opt,'skipIe')) opt.skipIe=false; end
if(~isfield(opt,'maxIeSteps')) opt.maxIeSteps=20; end
if(~isfield(opt,'etaDifCost')) opt.etaDifCost=1e-2; end
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

y=y(:);
Imea=exp(-y); alpha=xInit(:);

if(isfield(opt,'trueAlpha'))
    switch opt.errorType
        case 0
            trueAlpha = opt.trueAlpha/pNorm(opt.trueAlpha);
            computError= @(xxx) 1-(innerProd(xxx,trueAlpha)^2)/sqrNorm(xxx);
            clear 'trueAlpha'
        case 1
            trueAlphaNorm=sqrNorm(opt.trueAlpha);
            if(trueAlphaNorm==0) trueAlphaNorm=eps; end
            computError = @(xxx) sqrNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
            clear 'trueAlphaNorm'
        case 2
            trueAlphaNorm=pNorm(opt.trueAlpha);
            if(trueAlphaNorm==0) trueAlphaNorm=eps; end
            computError = @(xxx) pNorm(xxx-opt.trueAlpha)/trueAlphaNorm;
            clear 'trueAlphaNorm'
    end
end

if(opt.debugLevel>=2)
    figLin=figure;
    if(isfield(opt,'trueAlpha')) PhiTrueAlpha=Phi(opt.trueAlpha); end
end
if(opt.debugLevel>=3) figCost=figure; end
if(opt.debugLevel>=4) figRes=figure; end
if(opt.debugLevel>=5 && ~opt.skipIe) figIe=figure; end
if(opt.debugLevel>=6 && ~opt.skipAlpha) figAlpha=figure; end


switch lower(opt.sampleMode)
    case 'assigned'
        Ie=zeros(length(opt.kappa),1);
        temp=opt.kappa;
        temp1=abs(temp-1);
        temp2=find(temp1==min(temp1));
        Ie(temp2-1:temp2+1)=1/3;
    case 'logspan'
        kappa=logspace(-floor(opt.E/2)/(opt.E-1)*opt.logspan,...
            floor(opt.E/2-0.5)/(opt.E-1)*opt.logspan,opt.E);
        opt.E=length(kappa);
        Ie=zeros(opt.E,1); Ie(floor(opt.E/2)+1)=1;
end

switch lower(opt.IeStep)
    case lower('ActiveSet')
        IeStep = ActiveSet(B,b,Ie,opt.maxIeSteps,opt.stepShrnk);
    case lower('NPG')
        IeStep = NPG_ind(Ie,true,[],[],opt.maxIeSteps,opt.stepShrnk,opt.thresh);
    case lower('LBFGSB')
        IeStep=LBFGSB(Ie,opt.maxIeSteps);
end

switch lower(opt.noiseType)
    case lower('Gaussian')
        alphaStepFunc = @(III,aaa,pIout) gaussLAlpha(Imea,III,aaa,Phi,Phit,pIout,IeStep);
        IeStepFunc = @(A,III) gaussLI(Imea,A,III);
    case lower('Poisson')
        alphaStepFunc = @(III,aaa,pIout) poissLAlpha(Imea,III,aaa,Phi,Phit,pIout);
        IeStepFunc = @(A,III) poissLI(Imea,A,III);
end

q=kappa(2)/kappa(1);
polymodel=Spline(opt.spectBasis,[kappa(1)/q; kappa(:); kappa(end)*q]);
polyIout = polymodel.polyIout;

if(isfield(opt,'kappa'))
    polymodel.setPlot(opt.kappa,opt.iota,opt.epsilon);
    [opt.upkappa,opt.upiota]=getUpiota(opt.epsilon,opt.kappa,opt.iota);
end

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
    case {lower('NPGs'), lower('NPG')}
        switch(lower(opt.proximal))
            case lower('wvltADMM')
                proxmalProj=@(x,u,innerThresh,maxInnerItr) NPG.ADMM(Psi,Psit,x,u,...
                    innerThresh,maxInnerItr,false);
                penalty = @(x) pNorm(Psit(x),1);
                fprintf('Use l1 norm of wavelet coeff, ADMM\n');
            case lower('wvltLagrangian')
                proxmalProj=@(x,u,innerThresh,maxInnerItr) constrainedl2l1denoise(...
                    x,Psi,Psit,u,0,1,maxInnerItr,2,innerThresh,false);
                penalty = @(x) pNorm(Psit(x),1);
                fprintf('Use l1 norm of wavelet coeff, SPIRAL\n');
            case lower('tvl1')
                proxmalProj=@(x,u,innerThresh,maxInnerItr) Utils.denoiseTV(x,u,...
                    innerThresh,maxInnerItr,opt.mask,'l1');
                penalty = @(x) tlv(x,'l1');
                fprintf('Use l1 TV\n');
            case lower('tviso')
                proxmalProj=@(x,u,innerThresh,maxInnerItr) Utils.denoiseTV(x,u,...
                    innerThresh,maxInnerItr,opt.mask,'iso');
                penalty = @(x) tlv(x,'iso');
                fprintf('Use ISO TV\n');
        end

        if(strcmpi(opt.alphaStep,'npgs'))
            alphaStep=NPGs(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
            alphaStep.fArray{3} = penalty;
        elseif(strcmpi(opt.alphaStep,'npg'))
            alphaStep=NPG (1,alpha,opt.maxAlphaSteps,opt.stepShrnk,proxmalProj);
            alphaStep.fArray{3} = penalty;
        end
end
alphaStep.fArray{1} = @(aaa) alphaStepFunc(Ie,aaa,polyIout);
alphaStep.fArray{2} = @Utils.nonnegPen;
alphaStep.coef(1:2) = [1; opt.nu;];

B=eye(opt.E); b=zeros(opt.E,1);
if(isfield(opt,'CenterB') && opt.CenterB)
    if(~isfield(opt,'correctCenterB')) opt.correctCenterB=true; end
    temp=-eye(opt.E); temp(floor(opt.E/2)+1,:)=[]; temp(:,floor(opt.E/2)+1)=1;
    temp=temp/sqrt(2);
    B=[B; temp]; b=[b; zeros(opt.E-1,1)];
end

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

disp(['Use initial sparsity regulator u: ' num2str(alphaStep.u)]);

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

if(opt.debugLevel>=1)
    fprintf('%s\n', repmat( '=', 1, 80 ) );
    str=sprintf('Beam Hardening Correction %s_%s',opt.alphaStep,opt.IeStep);
    fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
    fprintf('%s\n', repmat('=',1,80));
    str=sprintf( ' %5s','itr');
    if(~opt.skipAlpha) 
        str=sprintf([str ' %12s'],'Obj');
        if(opt.continuation || strcmpi(opt.uMode,'relative'))
            str=sprintf([str ' %6s'],'u');
        end
        if(isfield(opt,'trueAlpha'))
            str=sprintf([str ' %12g'], 'error');
        end
        str=sprintf([str ' %12s %4s'], 'difα', 'αSrh');
        str=sprintf([str ' %12s'], 'difOjbα');
        if(strcmpi(opt.noiseType,'Gaussian'))
            str=sprintf([str ' %23s'], 'zmf(min,max)');
        end
    end
    if(~opt.skipIe)
        str=sprintf([str ' %4s'],'Istp');
        str=sprintf([str ' %12s'],'difI');
        str=sprintf([str ' %12s'],'difObjI');
    end
    fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
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

    p=p+1;
    str=sprintf(' %5d',p);
    
    % start optimize over alpha
    if(~opt.skipAlpha) 
        alphaStep.fArray{1} = @(aaa) alphaStepFunc(Ie,aaa,polyIout);
        preCost=alphaStep.cost;
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
        out.difCost(p)=abs(out.cost(p)-preCost)/out.cost(p);
        alpha = alphaStep.alpha;

        str=sprintf([str ' %12g'],out.cost(p));
        if(opt.continuation || strcmpi(opt.uMode,'relative'))
            out.uRecord(p,:)=[opt.u,alphaStep.u];
            str=sprintf([str ' %6g'],alphaStep.u);
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
            str=sprintf([str ' %12g'], out.RMSE(p));
        end

        str=sprintf([str ' %12g %4d'], out.difAlpha(p), alphaStep.ppp);
        str=sprintf([str ' %12g'], out.difCost(p));
        if(strcmpi(opt.noiseType,'Gaussian'))
            str=sprintf([str ' (%10g,%10g)'], IeStep.zmf(1), IeStep.zmf(2));
            out.zmf(p,:)=IeStep.zmf(:)';
        end
    end
    % end optimizing over alpha
     
    % start optimizing over Ie
    if(~opt.skipIe)
        % update the object fuction w.r.t. Ie
        if(strcmpi(opt.noiseType,'Gaussian') && strcmpi(opt.IeStep,'NPG'))
            if(max(IeStep.zmf(:))>=1 || p<20)  IeStep.maxItr=opt.maxIeSteps*10;
            else IeStep.maxItr=opt.maxIeSteps; end
        end

        A = polyIout(Phi(alpha),[]);
        IeStep.func = @(Ie) IeStepFunc(A,Ie);
        IeStep.thresh=opt.etaDifCost*out.difCost(p);
        preCost=IeStep.cost;
        [~,IeStepCnt]=IeStep.main();

        alphaStep.cost=IeStep.cost+opt.u*alphaStep.fVal(3);
        out.llI(p) = IeStep.cost;

        switch lower(opt.IeStep)
            case lower('ActiveSet')
                out.IeSteps(p)= IeStep.stepNum;
                out.course{p} = IeStep.course;
            case lower('NPG')
                out.IeSearch(p) = IeStep.ppp;
                out.IeSteps(p) = IeStepCnt;
            case lower('LBFGSB')
                out.IeSteps(p) = IeStepCnt;
        end
        str=sprintf([str ' %4d'],out.IeSteps(p));

        out.difIe(p)=norm(IeStep.Ie(:)-Ie(:))/norm(Ie);
        out.difllI(p)=abs(out.llI(p)-preCost)/out.llI(p);
        Ie = IeStep.Ie;
        str=sprintf([str ' %12g'], out.difIe(p));
        if(p>1)
            str=sprintf([str ' %12g'],out.difllI(p));
        else
            str=sprintf([str ' %12s'],' ');
        end
    end
    % end optimizing over Ie
    
    if(opt.debugLevel>1)
        if(opt.debugLevel>=2)
            set(0,'CurrentFigure',figLin);
            idx = randi(length(A),1000,1);
            s = linspace(min(A),max(A),1000);
            plot(A(idx),-log(Imea(idx)),'g.');
            hold on;
            if(isfield(opt,'trueAlpha'))
                PhiTrueAlpha=PhiTrueAlpha*innerProd(A,PhiTrueAlpha)/sqrNorm(PhiTrueAlpha);
                plot(PhiTrueAlpha(idx),-log(Imea(idx)),'.');
            end
            plot(s,-log(polyIout(s,IeStep.Ie)),'r-');
            hold off;
        end
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
if(isfield(opt,'trueAlpha')) fprintf(', RMSE=%g',out.RMSE(end)); end
fprintf('\n');

if(opt.estIe)
end

end

