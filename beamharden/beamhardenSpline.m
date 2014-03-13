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
%   $Revision: 0.3 $ $Date: Wed 12 Mar 2014 09:02:06 PM CDT
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

tic;
interiorPointAlpha=0; prpCGAlpha=1;
interiorPointIe=0; activeSetIe=1;
if(~isfield(opt,'K')) opt.K=2; end
if(~isfield(opt,'E')) opt.E=5; end
if(~isfield(opt,'u')) opt.u=1e-4; end
if(~isfield(opt,'skipAlpha')) opt.skipAlpha=0; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
if(~isfield(opt,'skipIe')) opt.skipIe=0; end
% The range for mass attenuation coeff is 1e-2 to 1e4 cm^2/g
if(~isfield(opt,'muRange')) opt.muRange=[1e-2; 1e4]; end
if(~isfield(opt,'sampleMode')) opt.sampleMode='exponential'; end
if(~isfield(opt,'maxAlphaSteps')) opt.maxAlphaSteps=1; end
if(~isfield(opt,'showImg')) opt.showImg=0; end
if(~isfield(opt,'visible')) opt.visible=1; end
if(~isfield(opt,'alphaStep')) opt.alphaStep='NCG_PR'; end
if(~isfield(opt,'continuation')) opt.continuation=false; end

Imea=exp(-y); alpha=xInit(:); Ie=zeros(opt.E,1);

if(isfield(opt,'trueAlpha'))
    trueAlpha = opt.trueAlpha/norm(opt.trueAlpha);
end

if(opt.showImg)
    figRes=1000; figure(figRes);
    if(~opt.skipAlpha) figAlpha=1001; figure(figAlpha); else figAlpha=0; end
    if(~opt.skipIe) figIe=1002; figure(figIe); else figIe=0; end
else
    figRes=0; figAlpha=0; figIe=0;
end
figAlpha=0;

switch lower(opt.sampleMode)
    case 'uniform'
        temp=linspace(opt.muRange(1),opt.muRange(2),opt.E);
        Ie(floor(opt.E/2)-1:floor(opt.E/2)+1)=1/3;
    case 'exponential'
        temp=logspace(log10(opt.muRange(1)),log10(opt.muRange(2)),opt.E);
        temp1=abs(temp-1);
        Ie(temp1==min(temp1))=1;
    case 'assigned'
        Ie=zeros(length(opt.mu),1);
        temp=opt.mu;
        temp1=abs(temp-1);
        temp2=find(temp1==min(temp1));
        Ie(temp2-1:temp2+1)=1/3;
    case 'logspan'
        temp=logspace(-floor((opt.E-1)/2)/(opt.E-1)*opt.logspan,...
            floor(opt.E/2)/(opt.E-1)*opt.logspan,opt.E);
        Ie(floor(opt.E/2+0.5))=1;
        if(strcmp(opt.spectBasis,'b0')) % extend to bigger end
            temp = [temp(:); temp(end)^2/temp(end-1)];
        elseif(strcmp(opt.spectBasis,'b1'))
            temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
        end
end

for i=1:opt.K-1
    mu(:,i)=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(opt.K-1)*eps);
end

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
        if(strcmp(opt.spectBasis,'dis'))
            Ie=interp1(opt.trueKappa, opt.trueUpiota,mu(:),'spline');
            temp = [mu(1);(mu(1:end-1)+mu(2:end))/2;mu(end)];
            temp = temp(2:end)-temp(1:end-1);
            Ie = Ie.*temp;
        elseif(strcmp(opt.spectBasis,'b0'))
            Ie=interp1(opt.trueKappa, opt.trueUpiota,mu(1:end-1),'spline');
        elseif(strcmp(opt.spectBasis,'b1'))
            Ie=interp1(opt.trueKappa, opt.trueUpiota,mu(2:end-1),'spline');
        end
        % there will be some points interplated negative and need to be removed
        Ie(Ie<0)=0;
    end
end

polymodel = Spline(opt.spectBasis,mu);
polymodel.setPlot(opt.trueKappa,opt.trueIota,opt.epsilon);
polyIout = polymodel.polyIout;

% find the best intial Ie starts
% R = polyIout(Phi(alpha),[]);
% for i=1:size(R,2)
%     temp(i) = var(y+log(R(:,i)),1);
% end
% idx = find(temp==min(temp));
% Ie = Ie*0;
% Ie(idx) = exp(-mean(y+log(R(:,idx))));

p=0; thresh=1e-4; str=''; strlen=0;
t1=0; thresh1=1e-8;
t2=0; thresh2Lim=1e-10;
if(interiorPointIe) 
    thresh2=1; t2Lim=1e-10;
else thresh2=1e-8; end

%max(Imea./(exp(-atten(Phi,alpha)*mu')*Ie))
if(interiorPointAlpha) 
    nonneg=@nonnegLogBarrier; 
else 
    nonneg=@nonnegPen; 
end
switch lower(opt.alphaStep)
    case lower('NCG_PR')
        alphaStep = Methods(3,alpha);
        alphaStep.fArray{2} = nonneg;
        if(isfield(opt,'muLustig'))
            fprintf('use lustig approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) lustigL1(aaa,opt.muLustig,Psi,Psit);
        end
        if(isfield(opt,'muHuber'))
            fprintf('use huber approximation for l1 norm\n');
            alphaStep.fArray{3} = @(aaa) huber(aaa,opt.muHuber,Psi,Psit);
        end
        alphaStep.coef(1:2) = [1; 1];
        alphaStep.maxStepNum = opt.maxAlphaSteps;
        alphaStep.stepShrnk = opt.stepShrnk;
        out.fVal =zeros(opt.maxItr,3);
    case {lower('SpaRSA')}
        alphaStep = Methods(2,alpha);
        alphaStep.coef(1:2) = [1; 1];
        alphaStep.fArray{2} = nonneg;
        alphaStep.maxStepNum = opt.maxAlphaSteps;
        alphaStep.stepShrnk = opt.stepShrnk;
        alphaStep.Psi = Psi;
        alphaStep.Psit = Psit;
        alphaStep.M = 5;
        out.fVal =zeros(opt.maxItr,3);
    case lower('FISTA_l1')
        alphaStep=FISTA_l1(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        alphaStep.coef(1:2) = [1; 1];
        alphaStep.fArray{2} = nonneg;
        out.fVal=zeros(opt.maxItr,3);
    case {lower('ADMM')}
        alphaStep=ADMM(1,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        out.fVal=zeros(opt.maxItr,2);
    case {lower('ADMM_N')}
        alphaStep = ADMM_N(2,alpha,opt.maxAlphaSteps,opt.stepShrnk,Psi,Psit);
        alphaStep.coef(1:2) = [1; 1];
        alphaStep.fArray{2} = nonneg;
        out.fVal=zeros(opt.maxItr,3);
end
out.llI         =zeros(opt.maxItr,1);
out.difAlpha    =zeros(opt.maxItr,1);
out.cost        =zeros(opt.maxItr,1);
out.course      = cell(opt.maxItr,1);
out.time        =zeros(opt.maxItr,1);
out.IeSteps     =zeros(opt.maxItr,1);
out.RMSE        =zeros(opt.maxItr,1);
out.difObjAlpha =zeros(opt.maxItr,1);
out.difObjIe    =zeros(opt.maxItr,1);

PsitPhitz=Psit(Phit(y));
PsitPhit1=Psit(Phit(ones(length(y),1)));

if(interiorPointIe)
    Ie(Ie<eps)=eps;
    while(sum(Ie)>1-eps)
        delta=sum(Ie)-(1-eps);
        temp=find(Ie>eps);
        numPos=length(temp);
        Ie(temp)=Ie(temp)-min( min(Ie(temp))-eps, delta/numPos  );
    end
else
    temp = polyIout(0,[]);
    B=[eye(opt.E); -temp(:)'/norm(temp)]; b=[zeros(opt.E,1); -1/norm(temp)];
    IeStep = ActiveSet(B,b,Ie);
    IeStep.maxStepNum = opt.maxIeSteps;
    IeStep.stepShrnk = opt.stepShrnk;
end

if(opt.continuation)
    [temp,temp1]=polyIout(0,Ie);
    t3=max(abs(PsitPhitz+PsitPhit1*log(temp)))*temp1/temp;
    alphaStep.u = t3*0.1;
else
    alphaStep.u = opt.u;
end

while( ~((alphaStep.converged || opt.skipAlpha) && (IeStep.converged || opt.skipIe)) )
    p=p+1;
    str=sprintf([str '\np=%-4d '],p);
    
    % start optimize over alpha
    if(~opt.skipAlpha)
        alphaStep.fArray{1} = @(aaa) gaussLAlpha(Imea,Ie,aaa,Phi,Phit,polyIout,IeStep);
        alphaStep.fArray{1} = @(aaa) linearModel(aaa,Phi,Phit,y);
        alphaStep.main();
        
        out.fVal(p,:) = (alphaStep.fVal(:))';
        out.cost(p) = alphaStep.cost;
        out.difAlpha(p) = norm(alphaStep.alpha(:)-alpha(:))^2;
        %out.alphaSearch(p) = alphaStep.ppp;
        alpha = alphaStep.alpha;
        if(opt.continuation && p>1 && ...
                abs(out.cost(p)-out.cost(p-1))/out.cost(p)<1e-3 && ...
                alphaStep.u>opt.u)
            alphaStep.u = max(alphaStep.u*0.2,opt.u);
            fprintf('\nnew u= %g\n',alphaStep.u);
            alphaStep.warned = true;
        end
        %if(out.stepSz~=s1) fprintf('lineSearch is useful!!\n'); end
        if(isfield(opt,'trueAlpha'))
            out.RMSE(p)=1-(alpha'*trueAlpha/norm(alpha))^2;
        end

        str=sprintf([str 'cost=%-10g RSE=%-10g ',...
            'dAlpha=%-10g aSearch=%d '],...
            out.cost(p),out.RMSE(p), out.difAlpha(p), ...
            alphaStep.ppp);
        if(p>1)
            str=sprintf([str 'pdObjAlpha=%g%% '],...
                (out.cost(p-1)-out.cost(p))/out.cost(p)*100);
        end
    end
    % end optimizing over alpha
    
    %if(out.delta<=1e-4) maxPP=5; end
    if(~opt.skipIe && ((~opt.skipAlpha && max(IeStep.zmf(:))<1) || (opt.skipAlpha)))
        % update the object fuction w.r.t. Ie
        IeStep.func = @(III) gaussLI(Imea,polyIout(Phi(alpha),[]),III);
        for i=1:opt.maxIeSteps
            IeStep.main();
        end
        Ie = IeStep.Ie;
        out.llI(p) = IeStep.cost;
        out.IeSteps(p)= IeStep.stepNum;
        out.course{p} = IeStep.course;
        %out.IeSearch(p) = IeStep.ppp;
        str=sprintf([str 'dObjIe=%-10g zmf=(%g,%g) IeSteps=%-3d'],...
            out.difObjIe(p),...
            IeStep.zmf(1), IeStep.zmf(2), out.IeSteps(p));
    end

    if(opt.showImg && p>1)
        set(0,'CurrentFigure',figRes);
        if(~opt.skipAlpha)
            subplot(2,1,1);
            semilogy(p-1:p,out.fVal(p-1:p,1),'r'); hold on;
            if(length(out.fVal(p,:))>1)
                semilogy(p-1:p,out.fVal(p-1:p,2),'g');
            end
            if(length(out.fVal(p,:))>2)
                semilogy(p-1:p,out.fVal(p-1:p,3),'b');
            end
        end
        if(~opt.skipIe)
            semilogy(p-1:p,out.llI(p-1:p),'c'); hold on;
            if(isfield(out,'penIe'))
                semilogy(p-1:p,out.penIe(p-1:p),'m:');
            end
        end
        semilogy(p-1:p,out.cost(p-1:p),'k');
        title(sprintf('cost(%d)=%g',p,out.cost(p)));

        if(~opt.skipAlpha && isfield(opt,'trueAlpha'))
            subplot(2,1,2);
            semilogy(p-1:p,out.RMSE(p-1:p)); hold on;
            title(sprintf('RMSE(%d)=%g',p,out.RMSE(p)));
        end
        drawnow;
    end
    
    if(figIe && ~opt.skipIe)
        set(0,'CurrentFigure',figIe);
        polymodel.plotSpectrum(Ie);
        title(sprintf('int upiota d kappa = %g',polyIout(0,Ie)));
        drawnow;
    end
    
    if(~opt.skipAlpha && figAlpha)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        %showImgMask(Qmask-Qmask1/2,opt.mask);
        %title(['size of Q=' num2str(length(Q))]);
        if(~isempty(IeStep.zmf))
            title(sprintf('zmf=(%g,%g)', IeStep.zmf(1), IeStep.zmf(2)))
        end
        drawnow;
    end
    %if(mod(p,100)==1 && p>100) save('snapshotFST.mat'); end
    if(opt.visible && p>1)
        if(alphaStep.warned || IeStep.warned)
            fprintf('%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
        strlen = length(str);
        str='';
    end
    if(p >= opt.maxItr) break; end
    out.time(p)=toc;
end

out.fVal(p+1:end,:) = [];
out.llI(p+1:end)=[]; out.time(p+1:end)=[]; out.RMSE(p+1:end)=[];
out.IeSteps(p+1:end)=[];
out.course(p+1:end) = []; out.difAlpha(p+1:end)=[];
out.difObjAlpha(p+1:end)=[]; out.difObjIe(p+1:end)=[];
out.Ie=Ie; out.mu=mu; out.alpha=alpha; out.cpuTime=toc; out.p=p;

out.opt = opt;

%if(activeSetIe && ~opt.skipIe) out.ASactive=ASactive; end
out.t2=t2; out.t1=t1;

fprintf('\n');

end

function [f,g,h] = nonnegLogBarrier(alpha)
    %if(any(alpha(:)<=0)) f=eps^-1; alpha(alpha<=0)=eps;
    %else f=-sum(log(alpha(:)));
    %end
    %if(nargout>1) g = -1./alpha; h=1./alpha.^2; end
    f=-sum(log(alpha(:)));
    if(nargout>1)
        g=-1./alpha;
        h=1./(alpha.^2);
    end
end

function [f,g,h] = nonnegPen(alpha)
    temp=(alpha<0);
    f=alpha(temp)'*alpha(temp);
    if(nargout>=2)
        g=zeros(size(alpha));
        g(temp)=2*alpha(temp);
        if(nargout>=3)
            h = @(x,opt) hessian(x,opt);
        end
    end
    function hh = hessian(x,opt)
        if(opt==1)
            hh = zeros(size(x));
            hh(temp,:) = x(temp,:)*2;
        else
            y = x(temp,:);
            hh = y'*y*2;
        end
    end
end

function [f,g,h]=barrierIe(Ie)
    %if(any(Ie)<=0)
    %    Ie(Ie<=0)=eps; f=eps^-1;
    %    if(1-sum(Ie)<=0) Ie=Ie*(1-eps)/sum(Ie); end
    %else
    %    if(1-sum(Ie)<=0) Ie=Ie*(1-eps)/sum(Ie); f=eps^-1;
    %    else f=-sum(log(Ie))-log(1-sum(Ie)); end
    %end
    f=-sum(log(Ie))-log(1-sum(Ie));
    if(nargout>1)
        g=1/(1-sum(Ie))-1./Ie;
        h=1/(1-sum(Ie))^2+diag(1./(Ie.^2));
    end
end

function [f,g,h] = lustigL1(alpha,xi,Psi,Psit)
    s=Psit(alpha);
    sqrtSSqrMu=sqrt(s.^2+xi);
    f=sum(sqrtSSqrMu);
    if(nargout>=2)
        g=Psi(s./sqrtSSqrMu);
        if(nargout>=3)
            h = @(x,opt) hessian(xi./(sqrtSSqrMu.^3),x,opt);
        end
    end
    function hh = hessian(weight,x,opt)
        y = Psit(x);
        if(opt==1)
            hh = Psi(weight.*y);
        else
            hh = y'*(weight.*y);
        end
    end
end

function [f,g,h] = huber(alpha,mu,Psi,Psit)
    s=Psit(alpha);
    idx = abs(s)<mu;
    temp = abs(s)-mu/2;
    temp(idx) = s(idx).^2/2/mu;
    f = sum(temp);
    if(nargout>=2)
        temp = ones(size(s));
        temp(idx) = s(idx)/mu;
        g=Psi(temp);
        if(nargout>=3)
            h = @(x,opt) hessian(x,opt);
        end
    end
    function hh = hessian(x,opt)
        y = Psit(x);
        if(opt==1)
            y(idx) = y(idx)/mu; y(~idx) = 0; hh=Psi(y);
        else
            y = y(idx); hh = y'*y/mu;
        end
    end
end

function [f,g,h] = linearModel(alpha,Phi,Phit,y)
    PhiAlpha=Phi(alpha);
    f=norm(y-PhiAlpha)^2/2;
    if(nargout>=2)
        g=Phit(PhiAlpha-y);
        if(nargout>=3)
            h=@(x,opt) hessian(x,opt);
        end
    end
    function hh = hessian(x,opt)
        yy = Phi(x);
        if(opt==1)
            hh = Phit(yy);
        else
            hh = yy'*yy;
        end
    end
end

