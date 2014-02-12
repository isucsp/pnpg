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
%   $Revision: 0.3 $ $Date: Tue 11 Feb 2014 09:42:00 PM CST
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
stepShrnk = 0.8;
interiorPointAlpha=0; prpCGAlpha=1;
interiorPointIe=0; activeSetIe=1;
if(~isfield(opt,'K')) opt.K=2; end
if(~isfield(opt,'E')) opt.E=5; end
if(~isfield(opt,'showImg')) opt.showImg=0; end
if(~isfield(opt,'skipAlpha')) opt.skipAlpha=0; end
if(~isfield(opt,'skipIe')) opt.skipIe=0; end
% The range for mass attenuation coeff is 1e-2 to 1e4 cm^2/g
if(~isfield(opt,'muRange')) opt.muRange=[1e-2; 1e4]; end
if(~isfield(opt,'sampleMode')) opt.sampleMode='exponential'; end
if(~isfield(opt,'visible')) opt.visible==1; end

Imea=exp(-y); alpha=xInit(:); Ie=zeros(opt.E,1);

if(isfield(opt,'trueAlpha'))
    opt.trueAlpha = opt.trueAlpha/norm(opt.trueAlpha);
end

if(opt.showImg)
    figRes=1000; figAlpha=1001; figIe=1002;
    figure(figAlpha); figure(figIe); figure(figRes);
else
    figRes=0; figAlpha=0; figIe=0;
end

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
        if(strcmp(opt.spectBasis,'b0'))
            temp = [temp(1)^2/temp(2); temp(:)];
        elseif(strcmp(opt.spectBasis,'b1'))
            temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
        end
end

for i=1:opt.K-1
    mu(:,i)=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(opt.K-1)*eps);
end

deltaEpsilon=mean([opt.epsilon(:) [opt.epsilon(2:end); opt.epsilon(end)]],2)-...
    mean([opt.epsilon(:) [opt.epsilon(1); opt.epsilon(1:end-1)]],2);
opt.trueIota=opt.trueIota/(opt.trueIota'*deltaEpsilon);

polymodel = Spline(opt.spectBasis,mu);
polymodel.setPlot(opt.trueKappa,opt.trueIota,opt.epsilon);
polyIout = polymodel.polyIout;

% find the best intial Ie starts
R = polyIout(Phi(alpha),[]);
for i=1:size(R,2)
    temp(i) = var(y+log(R(:,i)),1);
end
idx = find(temp==min(temp));
Ie = Ie*0;
Ie(idx) = exp(-mean(y+log(R(:,idx))));

% find the best intial Ie ends
if(opt.skipIe)
    Ie=interp1(opt.trueKappa(1:end-1),...
        abs(opt.trueIota(1:end-1)...
        .*(opt.epsilon(2:end)-opt.epsilon(1:end-1))...
        ./(opt.trueKappa(2:end)-opt.trueKappa(1:end-1))), ...
        mu(:),'spline');
    Ie(Ie<0)=0;
end
if(isfield(opt,'Ie') && length(opt.Ie)==length(mu(:))) Ie=opt.Ie(:); end;

p=0; thresh=1e-4; str='';
t1=0; thresh1=1e-8;
t2=0; thresh2Lim=1e-10;
if(interiorPointIe) 
    thresh2=1; t2Lim=1e-10;
else thresh2=1e-8; end

out.llAlpha=zeros(opt.maxItr,1);
out.nonneg=zeros(opt.maxItr,1);
out.llI=zeros(opt.maxItr,1);
out.cost=zeros(opt.maxItr,1);
out.course = cell(opt.maxItr,1);
out.time=zeros(opt.maxItr,1);
out.IeSteps = zeros(opt.maxItr,1);
out.RMSE=zeros(opt.maxItr,1);
out.deltaNormAlpha=zeros(opt.maxItr,1);
out.deltaNormIe=zeros(opt.maxItr,1);
out.llAlphaDif=zeros(opt.maxItr,1);

%max(Imea./(exp(-atten(Phi,alpha)*mu')*Ie))
llAlpha = @(aaa,III) gaussLAlpha(Imea,III,aaa,mu,Phi,Phit,polyIout);
llI = @(AAA,III) gaussLI(Imea,AAA,III);

if(interiorPointAlpha) 
    nonneg=@nonnegLogBarrier; 
else 
    nonneg=@nonnegPen; 
end
if(prpCGAlpha)
    alphaStep = ConjGrad(3);
    alphaStep.fArray{2} = nonneg;
    alphaStep.fArray{3} = @lustigL1;
    alphaStep.coef(1:2) = [1; 1];
end
if(isfield(opt,'a'))
    PsitPhitz=Psit(Phit(y));
    PsitPhit1=Psit(Phit(ones(length(y),1)));
end

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
    if(B(end,:)*Ie<b(end)) Ie=b(end)/(B(end,:)*Ie)*Ie; end
    Q = (B*Ie-b<1e-14);
    Z = null(B(Q,:),'r');

    IeStep = ActiveSet(@(III) llI(polyIout(Phi(alpha)),III), B,b,Ie);
    IeStep.maxStepNum = opt.maxIeSteps;
end

while( ~((alphaStep.converged || opt.skipAlpha) && (IeStep.converged || opt.skipIe)) )
    p=p+1;
    
    % start optimize over alpha
    if(~opt.skipAlpha)
        if(~isfield(opt,'t3'))
            [temp,temp1]=polyIout(0,Ie);
            t3=max(abs(PsitPhitz+PsitPhit1*log(temp)))*temp1/temp;
            t3 = t3*10^opt.a;
        end
        alphaStep.fArray{1} = @(xxx) llAlpha(xxx,Ie);
        alphaStep.coef(3) = t3;
        alphaStep.prCG();
        
        alpha = alphaStep.alpha;
        out.deltaNormAlpha(p)=deltaNormAlpha;
        out.t3(p) = t3;
        
        %if(out.stepSz~=s1) fprintf('lineSearch is useful!!\n'); end
        if(isfield(opt,'trueAlpha'))
            out.RMSE(p)=1-(alpha'*opt.trueAlpha/norm(alpha))^2;
        end
    end
    % end optimizing over alpha
    
    %if(out.delta<=1e-4) maxPP=5; end
    if(((~opt.skipAlpha && max(zmf(:))<1) || (opt.skipAlpha)) && ~opt.skipIe)
        % update the object fuction w.r.t. Ie
        IeStep.func = @(III) llI(polyIout(Phi(alpha),[]),III);
        IeStep.main();
        Ie = IeStep.Ie;
        out.llI(p) = IeStep.cost;
        out.IeSteps(p)= IeStep.stepNum;
        out.course{p} = IeStep.course;
        out.deltaNormIe(p) = IeStep.deltaNormIe;
    end
    
    if(out.llI(p)~=0) out.cost(p) = out.llI(p);
    else out.cost(p) = out.llAlpha(p); end
    if(~opt.skipAlpha && (isfield(out,'nonneg')))
        out.cost(p) = out.cost(p) + out.nonneg(p);
    end
    if(~opt.skipIe && isfield(out,'penIe'))
        out.cost(p) = out.cost(p) + out.penIe(p);
    end
    if(opt.showImg && p>1)
        cost = 0;
        set(0,'CurrentFigure',figRes);
        if(~opt.skipAlpha)
            subplot(2,1,1);
            semilogy(p-1:p,out.llAlpha(p-1:p),'g'); hold on;
            if (isfield(out,'nonneg'))
                semilogy(p-1:p,out.nonneg(p-1:p),'b--');
            end
        end
        if(~opt.skipIe)
            semilogy(p-1:p,out.llI(p-1:p),'r'); hold on;
            cost = cost + out.llI(p);
            if(isfield(out,'penIe'))
                semilogy(p,out.penIe(p),'m--');
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
    
    if(figIe)
        set(0,'CurrentFigure',figIe);
        polymodel.plotSpectrum(Ie);
        title(sprintf('int upiota d kappa = %g',polyIout(0,Ie)));
        drawnow;
    end
    
    if(~opt.skipAlpha && figAlpha)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        %showImgMask(Qmask-Qmask1/2,opt.mask);
        %title(['size of Q=' num2str(length(Q))]);
        title(sprintf('zmf=(%g,%g)', zmf(1), zmf(2)))
        drawnow;
    end
    %if(mod(p,100)==1 && p>100) save('snapshotFST.mat'); end
    if(opt.visible)
        strlen = length(str);
        str=sprintf('\np=%-4d cost=%-10g RSE=%-10g dAlpha=%-10g dIe=%-10g zmf=(%g,%g) IeSteps=%-3d',...
            p,out.cost(p),out.RMSE(p), out.deltaNormAlpha(p), ...
            out.deltaNormIe(p), zmf(1),zmf(2), out.IeSteps(p));
        if(alphaStep.warned || IeSteps.warned)
            fprintf('%s',str);
        else
            fprintf([repmat('\b',1,strlen) '%s'],str);
        end
    end
    if(p >= opt.maxItr) break; end
    out.time(p)=toc;
end

out.llAlpha(p+1:end) = []; out.nonneg(p+1:end) = [];
out.llI(p+1:end)=[]; out.time(p+1:end)=[]; out.RMSE(p+1:end)=[];
out.llAlphaDif(p+1:end)=[]; out.IeSteps(p+1:end)=[];
out.course(p+1:end) = [];
out.deltaNormAlpha(p+1:end)=[]; out.deltaNormIe(p+1:end)=[];
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
            h = @(xl,xr) hessian(xl,xr);
        end
    end
    function hh = hessian(xl,xr)
        keyboard
        if(isempty(xl))
            y = xr(temp,:);
            hh = y'*y*2;
        elseif(xl==1)
            hh = zeros(size(xr));
            hh(temp,:) = xr(temp,:)*2;
        else
            hh = xl(temp,:)'*xr(temp,:)*2;
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

function [f,g,h] = lustigL1(alpha,xi,Psit)
    s=Psit(alpha);
    sqrtSSqrMu=sqrt(s.^2+xi);
    f=sum(sqrtSSqrMu);
    g=Psi(s./sqrtSSqrMu);
    h = @(x,opt) hessian(xi./(sqrtSSqrMu.^3),x,opt);
    function hh = hessian(weight,x,opt)
        keyboard
        y = Psit(x);
        if(opt==1)
            hh = Psi(weight.*y);
        else
            hh = y'*(weight.*y);
        end
    end
end

