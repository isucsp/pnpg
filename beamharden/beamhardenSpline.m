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
%   $Revision: 0.3 $ $Date: Thu 30 Jan 2014 10:41:38 PM CST
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
%

tic;
K=2; E=5;
skipAlpha=0;
interiorPointAlpha=0;
prpCGAlpha=1;

if(isfield(opt,'skipIe')) skipIe=opt.skipIe; else skipIe=0; end
interiorPointIe=0;
activeSetIe=1;

if(isfield(opt,'K')) K=opt.K; end
if(isfield(opt,'E')) E=opt.E; end

Imea=exp(-y); alpha=xInit(:); Ie=zeros(E,1);

% The range for mass attenuation coeff is 1e-2 to 1e4 cm^2/g
if(isfield(opt,'muRange')) temp=opt.muRange(:); else temp=[1e-2; 1e4]; end
if(~isfield(opt,'sampleMode')) opt.sampleMode='exponential'; end
if(isfield(opt,'showImg') && opt.showImg==1) show=1; else show=0; end
if(isfield(opt,'visible') && opt.visible==1) visible=1; else visible=0; end
if(show)
    figRes=1000; figAlpha=1001; figIe=1002;
    figure(figAlpha); figure(figIe); figure(figRes);
else
    figRes=0; figAlpha=0; figIe=0;
end

switch lower(opt.sampleMode)
    case 'uniform'
        temp=linspace(temp(1),temp(2),E);
        Ie(floor(E/2)-1:floor(E/2)+1)=1/3;
    case 'exponential'
        temp=logspace(log10(temp(1)),log10(temp(2)),E);
        temp1=abs(temp-1);
        Ie(temp1==min(temp1))=1;
    case 'assigned'
        Ie=zeros(length(opt.mu),1);
        temp=opt.mu;
        temp1=abs(temp-1);
        temp2=find(temp1==min(temp1));
        Ie(temp2-1:temp2+1)=1/3;
    case 'logspan'
        temp=logspace(-floor((E-1)/2)/(E-1)*opt.logspan,...
            floor(E/2)/(E-1)*opt.logspan,E);
        Ie(floor(E/2+0.5))=1;
        if(strcmp(opt.spectBasis,'b0'))
            temp = [temp(1)^2/temp(2); temp(:)];
        elseif(strcmp(opt.spectBasis,'b1'))
            temp = [temp(1)^2/temp(2); temp(:); temp(end)^2/temp(end-1)];
        end
end

for i=1:K-1
    mu(:,i)=temp(:);  %*mean(X(find(idx(:)==i+1))); %/(1-(K-1)*eps);
end

switch lower(opt.spectBasis)
    % for b0-spline, length(I)=length(kappa)-1;
    case 'b0'       % B-0 spline with nodes be kappa
        polyIout = @b0Iout;
        plotSpectrum = @(III) plotB0Upiota(opt.trueMu(1:end-1),...
            abs(opt.trueIe(1:end-1)...
            .*(opt.epsilon(2:end)-opt.epsilon(1:end-1))...
            ./(opt.trueMu(2:end)-opt.trueMu(1:end-1))), ...
            mu, III);
    case 'dis'
        polyIout = @disIout;
        plotSpectrum = @(III) plotDisUpiota(opt.trueMu(1:end-1),...
            abs(opt.trueIe(1:end-1)...
            .*(opt.epsilon(2:end)-opt.epsilon(1:end-1))...
            ./(opt.trueMu(2:end)-opt.trueMu(1:end-1))), ...
            mu, III);
end

% find the best intial Ie starts
R = (polyIout(mu,Phi(alpha)));
for i=1:size(R,2)
    temp(i) = var(y+log(R(:,i)),1);
end
idx = find(temp==min(temp));
Ie = Ie*0;
Ie(idx) = exp(-mean(y+log(R(:,idx))));

% find the best intial Ie ends
if(skipIe)
    Ie=interp1(opt.trueMu,opt.trueIe,mu(:),'spline');
    Ie(Ie<0)=0;
    epsilon=interp1(opt.trueMu,opt.epsilon,mu(:),'spline');
    deltaEpsilon=mean([epsilon(:) [epsilon(2:end); epsilon(end)]],2)-...
        mean([epsilon(:) [epsilon(1); epsilon(1:end-1)]],2);
    Ie=Ie.*abs(deltaEpsilon);
end
if(isfield(opt,'Ie') && length(opt.Ie)==length(mu(:))) Ie=opt.Ie(:); end;
if(isfield(opt,'skipAlpha') && opt.skipAlpha==1) skipAlpha=1; end;
if(isfield(opt,'t3'))
    t3=opt.t3; %abs(costA/costLustig)*1e-3;
    out.t3=t3;
end
if(isfield(opt,'a'))
    PsitPhitz=Psit(Phit(y));
    PsitPhit1=Psit(Phit(ones(length(y),1)));
end

if(prpCGAlpha) preP=0; preG=1; end
if(activeSetIe) 
    minZHZ=0; end

alphaReady=0; IeReady=0;
p=0; maxItr=opt.maxItr; thresh=1e-4; str=[];
t1=0; thresh1=1e-8;
t2=0; thresh2Lim=1e-10;
stepShrnk=0.9;
if(interiorPointIe) 
    thresh2=1; t2Lim=1e-10; else thresh2=1e-8; end

out.llAlpha=zeros(maxItr,1);
out.penAlpha=zeros(maxItr,1);
out.llI=zeros(maxItr,1);
out.time=zeros(maxItr,1);
out.IeSteps = zeros(maxItr,1);
out.RMSE=zeros(maxItr,1);
out.llAlphaDif=zeros(maxItr,1);
%ASactive=zeros(maxItr,1);
asIdx=0;

muLustig=opt.muLustig;

%max(Imea./(exp(-atten(Phi,alpha)*mu')*Ie))
llAlpha = @(aaa,III) gaussLAlpha(Imea,III,aaa,mu,Phi,Phit,polyIout);
llI = @(AAA,III) gaussLI(Imea,AAA,III);

if(interiorPointAlpha) 
    penAlpha=@barrierAlpha; 
else 
    penAlpha=@barrierAlpha2; 
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
    temp = polyIout(mu,0);
    B=[eye(E); -temp(:)'/norm(temp)]; b=[zeros(E,1); -1/norm(temp)];
    if(B(end,:)*Ie<b(end)) Ie=b(end)/(B(end,:)*Ie)*Ie; end
    Q = (B*Ie-b<1e-14);
    Z = null(B(Q,:),'r');

    IeStep = ActiveSet(B,b,Ie);
end

while( ~((alphaReady || skipAlpha) && (IeReady || skipIe)) )
    p=p+1;
    
    % start optimize over alpha
    if(~skipAlpha)
        [costA,zmf,diff0,weight]=llAlpha(alpha,Ie);
        [costB,difphi,hphi]=penAlpha(alpha);
        
        if(min(weight)<=0)
            fprintf(['\nfAlpha[warning]: obj is non-convex over alpha,'...
                'minimum=%g\n'],min(weight));
            str='';
        end
        
        s=Psit(alpha);
        sqrtSSqrMu=sqrt(s.^2+muLustig);
        costLustig=sum(sqrtSSqrMu);
        difLustig=Psi(s./sqrtSSqrMu);
        
        if(t1==0 && p==1)
            if(interiorPointAlpha)
                t1=min([1, abs(diff0'*difphi/norm(difphi)), abs(costA/costB)]);
            else t1=1;
            end
        end
        
        if(~isfield(opt,'t3'))
            [temp,temp1]=polyIout(mu,0,Ie);
            t3=max(abs(PsitPhitz-PsitPhit1*log(sum(Ie))))*temp1/temp;
            t3=t3*10^opt.a; out.t3 = t3;
        end
        cost=costA+t1*costB+t3*costLustig;
        difAlpha=diff0+t1*difphi+t3*difLustig;
        if(0)
            %afun=@(xxx,yyy) fhessianA(xxx);
            %[deltaAlpha,~]=bicg(afun,difAlpha,[],5);
            %afun=@(xxx,yyy) fhessianA(xxx);
            fhessianA=@(gAlpha) hessianA(gAlpha,weight,hphi*t1,Phi,Phit);
            fatHessianA=@(gAlpha) atHessianA(gAlpha,weight,t1*hphi,Phi,Phit);
            deltaAlpha = cg(-difAlpha,fhessianA,fatHessianA,1);
        end
        if(prpCGAlpha)
            beta=difAlpha'*(difAlpha-preG)/(preG'*preG);
            deltaAlpha=difAlpha+max(beta,0)*preP;
            normDelta=difAlpha'*deltaAlpha;
            s1=normDelta/atHessianA(deltaAlpha,weight,t1*hphi,Phi,Phit,...
                t3, Psit,muLustig,sqrtSSqrMu);
            preP=deltaAlpha; preG=difAlpha;
            deltaAlpha=deltaAlpha*s1;
        end
        
        if(interiorPointAlpha)
            temp=find(deltaAlpha>0);
            if(isempty(temp)) maxStep=1;
            else maxStep=min(alpha(temp)./deltaAlpha(temp)); end
            maxStep=maxStep*0.99;
        else maxStep=1;
        end
        
        penalty=@(x) t1*penAlpha(x)+t3*sum(sqrt(Psit(x).^2+muLustig));
        
        % start of line Search
        pp=0; stepSz=min(1,maxStep);
        while(1)
            pp=pp+1;
            newX=alpha-stepSz*deltaAlpha;
            %newX(newX<0)=0; % force it be positive;
            
            [newCostA,zmf]=llAlpha(newX,Ie);
            [newCostB]=penalty(newX);
            newCost=newCostA+newCostB;
            
            if(newCost <= cost - stepSz/2*difAlpha'*deltaAlpha)
                break;
            else
                if(pp>10) stepSz=0; else stepSz=stepSz*stepShrnk; end
            end
        end
        %end of line search
        
        out.llAlphaDif(p) = norm(alpha(:)-newX(:))^2;
        out.llAlpha(p)=newCostA; out.penAlpha(p) = newCostB;
        out.time(p)=toc;
        
        alpha = newX;
        deltaNormAlpha=difAlpha'*deltaAlpha;
        
        %if(out.stepSz~=s1) fprintf('lineSearch is useful!!\n'); end
        if(interiorPointAlpha)
            if(deltaNormAlpha<1e-5)
                if(t1 < 1e-10/length(alpha))
                    alphaReady=1;
                else
                    t1=t1/10;
                    thresh1=thresh1/10;
                end
            end
        else
            if(deltaNormAlpha< thresh1)
                if(t1 < 1e2) 
                    t1=t1*10;
                else alphaReady=1;
                end
            else
                if(stepSz==0)
                    t1=max(1,t1/10);
                end
            end
        end
    end
    % end optimizing over alpha
    
    pp=0; maxPP=opt.maxIeStep; IeReady=false;
    %if(out.delta<=1e-4) maxPP=5; end
    A = polyIout(mu,Phi(alpha));
    %IeStep.main(A,Ie,maxPP);
    while(((~skipAlpha && max(zmf(:))<1) || (skipAlpha)) && pp<maxPP && ~skipIe)
        pp=pp+1;
        [costA,zmf,diff0,h0] = llI(A,Ie);
        if(interiorPointIe)
            optIeInteriorPoint; else optIeActiveSet; end
    end
    out.IeSteps(p)=pp;
    
    if(p >= maxItr) 
        break; end
    if(show)
        %figure(911); showImgMask(alpha,opt.mask);
        %figure; showImgMask(deltaAlpha,opt.mask);
        cost = 0;
        if(figRes)
            nCurve=0;
            set(0,'CurrentFigure',figRes);
            if(~skipAlpha && (isfield(out,'penAlpha')))
                semilogy(p,out.penAlpha(p),'b.');
                nCurve=nCurve+1;
                strLegend{nCurve}='alpha penalty';
                cost = cost + out.penAlpha(p);
            end
        end
        if(~skipIe)
            semilogy(p,out.llI(p),'r.'); hold on;
            nCurve = nCurve + 1;
            strLegend{nCurve}='likelihood';
            cost = cost + out.llI(p);
            if(isfield(out,'penIe'))
                semilogy(p,out.penIe(p),'m.');
                cost = cost + out.penIe(p);
                nCurve=nCurve+1;
                strLegend{nCurve}='Ie penalty';
            end
        end
        semilogy(p,cost,'k.');
        nCurve=nCurve+1;
        strLegend{nCurve}='overall cost';
        title(['total cost=' num2str(cost)]);
        legend(strLegend); drawnow;
    end
    
    if(figIe)
        set(0,'CurrentFigure',figIe);
        plotSpectrum(Ie);
        title(sprintf('int upiota d kappa = %g',polyIout(mu,0,Ie)));
        drawnow;
    end
    
    if(~skipAlpha && figAlpha)
        set(0,'CurrentFigure',figAlpha); showImgMask(alpha,opt.mask);
        %showImgMask(Qmask-Qmask1/2,opt.mask);
        %title(['size of Q=' num2str(length(Q))]);
        title(['zmf=' num2str(max(zmf))])
        drawnow;
    end
    if(0)
        fprintf(repmat('\b',1,length(str)));
        str=sprintf('\n%-8s%-13s%-13s%-13s%-13s%-13s%-13s\n',...
            'itr','normDif1','normDif2','t1','t2','s1','s2');
        str=[str, sprintf('%d-0:\t%-13g%-13g%-13g%-13g%-13g%-13g\n',...
            p,out.delta,out.delta,t1,t2,out.stepSz*s1,out.stepSz)];
        str=[str, sprintf('%-13s%-13s%13s-%-13s%-13s%-13s\n',...
            'cost','f0','(z','f)','penAlpha','f2')];
        str=[str, sprintf('%-13g%-13g%13g~%-13g%-13g%-13g\n',...
            out.cost, out.costA, out.zmf(1),...
            out.zmf(2), out.costB, out.costB)];
        fprintf('%s',str);
    end
    if(0)
        fprintf(repmat('\b',1,length(str)));
        str=sprintf('%d  %-13g%-13g',...
            p,t1,t2);
        str=[str, sprintf('%-13g%-13g%-13g',...
            out.cost, out.costA, ...
            out.zmf(2))];
        fprintf('%s',str);
    end
    %if(mod(p,100)==1 && p>100) save('snapshotFST.mat'); end
    out.RMSE(p)=1-(alpha'*opt.trueAlpha/norm(alpha))^2;
end
out.llAlpha(p+1:end) = []; out.penAlpha(p+1:end) = [];
out.llI(p+1:end)=[]; out.time(p+1:end)=[]; out.RMSE(p+1:end)=[];
out.llAlphaDif(p+1:end)=[]; out.IeSteps(p+1:end)=[];
out.Ie=Ie; out.mu=mu; out.alpha=alpha; out.cpuTime=toc; out.p=p;

if(activeSetIe && ~skipIe) 
    out.ASactive=ASactive; end
out.t2=t2; out.t1=t1;

fprintf('\n');

end

function [f,g,h] = barrierAlpha(alpha)
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

function [f,g,h] = barrierAlpha2(alpha)
    temp=(alpha<0);
    f=alpha(temp)'*alpha(temp);
    g=2*alpha;
    g(~temp)=0;
    h=2*ones(size(g)); h(~temp)=0;
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

function  h=atHessianA(gAlpha,weight,t1hphi,Phi,Phit,t3,Psit,muLustig,s1)
    temp=Phi(gAlpha);
    h=2*temp'*(weight.*temp);
    h=h+gAlpha'*(t1hphi.*gAlpha);
    if(nargin>5)
        temp=Psit(gAlpha);
        h=h+t3*muLustig*(s1.^(-3))'*(temp.^2);
    end
end

function h=hessianA(gAlpha,weight,t1hphi,Phi,Phit)
    temp=Phi(gAlpha);
    h=2*Phit(weight.*temp);
    h=h+(t1hphi.*gAlpha);
end

function x= cg(c,hessianA,atHessianA,maxItr)
    % This function solve the problem 
    % min c'*x+1/2 atHessianA(x)
    % hessianA=hessian*x
    x=0; g=c; p=0; i=0;
    while(i<maxItr)
        i= i+1;
        preP=p; preG=g;
        g=c+hessianA(x);
        p=-g-g'*g/(preG'*preG)*preP;
        x=x-p*(p'*g/atHessianA(p));
    end
end



function plotB0Upiota(trueMu, trueUpiota, mu, Ie)
    loglog(trueMu,trueUpiota,'r.-'); hold on;
    mu=reshape([mu(:)';mu(:)'],[],1);
    mu(1)=[]; mu(end)=[];
    Ie=reshape([Ie(:)';Ie(:)'],[],1);
    loglog(mu,Ie,'*-'); hold off;
    %ylim([1e-10 1]);
    xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
end

function plotDisUpiota(trueMu, trueUpiota, mu, Ie)
    loglog(trueMu,trueUpiota,'r.-'); hold on;
    temp=[mu(:); mu(end)^2/mu(end-1)];
    mu=reshape([temp';temp'],[],1);
    mu(1)=[]; mu(end)=[];
    temp=temp(2:end)-temp(1:end-1);
    Ie = Ie./temp;
    Ie=reshape([Ie(:)';Ie(:)'],[],1);
    loglog(mu,Ie,'*-'); hold off;
    %ylim([1e-10 1]);
    xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
end

