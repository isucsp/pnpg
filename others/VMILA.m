function out = ...
    VMILA(gn, H, HT, bg, rho, grad1, grad2, div, NIT,init, verbose,...
        obj,eta,P,p,alpha_min,alpha_max,inn_ini,thresh)

% ***************** * VMILA package * *****************
% 
% This software implements the Variable Metric Inexact Line?search
% Algorithm (VMILA) described in the paper 
% 
% S. Bonettini, I. Loris, F. Porta, M. Prato, Variable metric inexact
% line-search based methods for nonsmooth optimization
% 
% currently submitted for publication and available at the website
% http://arxiv.org/abs/1506.00385
% 
% The VMILA package includes four .m and three .mat files, which allow to
% reproduce the VMILA plots in the numerical experiments described in
% Section 6.1 of the paper.
% 
% In particular:
% 
% 1. main.m is the script with the choice of the dataset and the call to
% VMILA
% 
% 2. VMILA.m is the function implementing VMILA
% 
% 3. FISTA.m is the function implementing the Fast Iterative
% Shrinkage-Thresholding algorithm by Beck and Teboulle (SIAM J Imagin Sci
% 2 (2009), 183-202), which is used as inner solver in VMILA
% 
% 4. dctshift.m is a function creating an array containing the first column
% of a blurring matrix when implementing reflexive boundary conditions (see
% Chapter 4 of P.C. Hansen, J.G. Nagy, D.P. O'Leary, Deblurring Images -
% Matrices, Spectra, and Filtering, SIAM, Philadelphia, 2006).
% 
% As for the .mat files,
% 
% 1. TP_cameraman_RBC is the dataset based on the 256x256 Matlab cameraman
% image
% 
% 2. TP_fantoccio_RBC is the dataset based on the 256x256 Matlab
% Shepp-Logan phantom
% 
% 3. TP_micro_RBC is the dataset based on a 128x128 confocal microscopy
% phantom (see also R.M. Willet and R.D. Nowak, IEEE Trans. Med. Imaging 22
% (2003), 332?350)
% 
% The software is a beta version and further optimizations and
% documentation are still in progress.
% 
% For comments and further information please send us an email at
% silvia.bonettini@unife.it.


nobj = length(obj);
for i = 1 : nobj
    normobj(i) = norm(obj{i}(:));
    err{i}     = zeros(NIT+1,1); %arrays to store errors per iteration
end

Primal     = zeros(NIT+1,1);
alpha_vec  = zeros(NIT+1,1);
nit_vec    = zeros(NIT+1,1);
TimeCost   = zeros(NIT+1,1);
hsigma_vec = zeros(NIT+1,1);
G_vec      = zeros(NIT+1,1);

n = length(gn);
zeroindex = gn <= 0;
nonzeroindex = ~zeroindex;
ONE  = ones(n);
%rapp = zeros(n);
rapp = zeros(size(gn));         % modified by Renliang Gu
%HTe = HT(ONE);
HTe = HT(ones(size(gn)));       % modified by Renliang Gu
if min(HTe(:))== max(HTe(:))
    HTe = HTe(1);
end

gamma = 1e-4; % linesearch parameters
beta = 0.4;

Malpha = 3;   % memory on alphaBB
tau = 0.5;
Valpha = 1e30 * ones(Malpha,1);

soglia_X_low = 1e-10;
soglia_X_upp = 1e+10; % threshold on scaling matrix
alpha = 1.3;          % initial alpha

% initial point
%x = max(gn,0);
%x = max(HT(gn),0);              % modified by Renliang Gu
x = init;              % modified by Renliang Gu

% gradient and objective function at current x
dx1     = grad1(x);
dx2     = grad2(x);
densq   = sqrt(dx1.^2 + dx2.^2);
TV      = sum(densq(:));
den     = H(x) + bg;
rapp(nonzeroindex) = gn(nonzeroindex)./ den(nonzeroindex);
KL = sum( gn(nonzeroindex).* log(rapp(nonzeroindex)) + den(nonzeroindex) - gn(nonzeroindex) );
KL = KL + sum(den(zeroindex));
fv = KL + rho*TV; % objective function

g_KL = HTe - HT(rapp); 
g = g_KL; % gradient
Primal(1) = fv;
alpha_vec(1) = alpha;
for i = 1 : nobj
    err{i}(1) = norm(obj{i}(:)-x(:))/normobj(i);
end

if verbose
    fprintf('\nInitial: Primal=%8.3e', Primal(1));
    for i=1:nobj
        fprintf(' err{%g} %g', i, err{i}(1));
    end
end

% start CPU clock
TimeCost(1)=0;
t0 = tic;

% step 1 - scaling matrix
X =  x./HTe; 
X(X < soglia_X_low) = soglia_X_low;
X(X > soglia_X_upp) = soglia_X_upp;
D = 1./X; 

for itr = 1:NIT
    Valpha(1:Malpha-1) = Valpha(2:Malpha);

    % step 2 - compute ytilde
    z = x - alpha*X.*g_KL;
    if itr == 1 || inn_ini == 0
       v10 = zeros(size(x)); v20 = v10; v30 = v10;
    else
       v10 = v1; v20 = v2; v30 = v3;
    end
    inNIT = 1500;
    a = 2.1;
    gamma_inn = 1/(9*max(X(:)))/alpha;
    [y,v1,v2,v3,hsigma,G,nit,dy1,dy2] = FISTA(v10,v20,v30,gamma_inn,rho,...
        alpha,X,D,grad1,grad2,div,z,rho*TV,g_KL,inNIT,eta,a,max(0,verbose-1));
    hsigma_vec(itr) = hsigma;
    G_vec(itr) = G;
    nit_vec(itr) = nit;
    
    % step 3 - compute d
    d = y - x;
    gd = hsigma;
    denold = den;
    dx1old = dx1; dx2old = dx2;
    lam = 1;
    Hd  = H(d);
    A1d = dy1 - dx1;
    A2d = dy2 - dx2;

    % step 4 - linesearch
    fcontinue = 1;
    
    fr = fv;

    while fcontinue
        xplus = x + lam*d;
        
        % update the objective function value
        dx1     = dx1old + lam*A1d; 
        dx2     = dx2old + lam*A2d;
        densq   = sqrt(dx1.^2 + dx2.^2);
        TV      = sum(densq(:));
        den     = denold + lam*Hd;
        rapp(nonzeroindex) = gn(nonzeroindex)./ den(nonzeroindex);
        KL = sum( gn(nonzeroindex).* log(rapp(nonzeroindex)) + den(nonzeroindex) - gn(nonzeroindex) );
        KL = KL + sum(den(zeroindex));
        fv = KL + rho*TV;   % objective function at xplus
        
        % Armijo condition
        if ( fv <= fr + gamma * lam * gd || lam<1e-12)
            difX=relativeDif(x,xplus);
            x = xplus; clear xplus;
            sk = lam*d; % difference between iterates            
            g_KL    = HTe - HT(rapp);
            gtemp   = g_KL ;
            yk = gtemp - g; % difference between gradients
            g = gtemp; clear gtemp;
            fcontinue = 0;
        else
            lam = lam * beta;
        end
    end
    
    for i=1:nobj
        err{i}(itr + 1) = norm(obj{i}(:)-x(:))/normobj(i);
    end
    if verbose
        fprintf('\n%4d): f(x)=%g  alpha %g nit %g', itr, ...
            fv, alpha, nit );
        for i=1:nobj
            fprintf(' err{%g} %g', i, err{i}(itr + 1));
        end
    end
    
    % compute steplength and scaling matrix for the next iteration
    soglia_X_upp = sqrt(1 + P/itr^p);
    soglia_X_low = 1/soglia_X_upp;
    
    X =  x./HTe;
    X(X < soglia_X_low) = soglia_X_low;
    X(X > soglia_X_upp) = soglia_X_upp;
    D = 1./X;

    sk2 = sk.*D; yk2 = yk.*X;

    bk = sum(dot(sk2,yk));  ck = sum(dot(yk2,sk));
    if (bk <= 0)
        alpha1 = alpha_max;
    else
        alpha1BB = sum(dot(sk2,sk2))/bk;
        alpha1 = min(alpha_max, max(alpha_min, alpha1BB));
    end
    if (ck <= 0)
        alpha2 = alpha_max;
    else
        alpha2BB = ck/sum(dot(yk2,yk2));
        alpha2 = min(alpha_max, max(alpha_min, alpha2BB));
    end

    Valpha(Malpha) = alpha2;
    
    if (alpha2/alpha1 < tau)
        alpha = min(Valpha);
        tau = tau*0.9;
    else
        alpha = alpha1;
        tau = tau*1.1;
    end

    Primal(itr + 1)   = fv;
    TimeCost(itr + 1) = toc(t0);
    alpha_vec(itr + 1) = alpha;

    if(difX<thresh) break; end  % add by Renliang Gu
    
end

Primal(itr+2:end) = [];
nit_vec(itr+1) = [];
hsigma_vec(itr+1) = [];
G_vec(itr+1) = [];
alpha_vec(itr+2:end) = [];
TimeCost(itr+2:end) = [];

out.x=x;
out.time=TimeCost;
out.cost=Primal;
out.alpha=alpha_vec;
out.RMSE=err;
out.innerItr=nit_vec;
out.hsigma=hsigma_vec;
out.G=G_vec;

if verbose
    fprintf('\n');
end




function [ybar,w1,w2,w3,hsigmaplus,G,iter,g1ybar,g2ybar,normgybar,varargout]...
        = FISTA(v1,v2,v3,...
        gamma,beta,alpha,X,D,grad1,grad2,div,z,f1k,g,NIT,eta,a,varargin)

psivec = zeros(NIT+1,1);
hsigmavec = zeros(NIT+1,1);
if nargin == 18
    verbose = varargin{1};
else
    verbose = 2;
end
hsigmafixedpart = -f1k -alpha/2*(X(:)'*(g(:).^2));
Psifixedpart    = hsigmafixedpart + (D(:)'*(z(:).^2))/2/alpha;
alphaX = alpha*X;
w1old = v1; w2old = v2; w3old = v3;
ATv = div(v1,v2)+v3; ATw = ATv;
y   = z - alphaX.*ATv;
gstarv = 0;
Psi = -(D(:)'*(y(:).^2))/2/alpha - gstarv + Psifixedpart;
psivec(1) = -Psi;
g1y = grad1(y); g2y = grad2(y); 

ybar   = max(y,0);
g1ybar = grad1(ybar); g2ybar = grad2(ybar); normgybar = sqrt(g1ybar.^2+g2ybar.^2);
f1ybar = beta*sum(normgybar(:));
hsigmaplus = (D(:)'*(ybar(:)-z(:)).^2)/2/alpha + f1ybar + hsigmafixedpart;
hsigmavec(1) = hsigmaplus;
if verbose
    fprintf('\n      0) hsigmaplus %g psi %g gamma %g',hsigmaplus,Psi,gamma);
end

for iter = 1:NIT
    w1 = v1 + gamma*g1y;
    w2 = v2 + gamma*g2y;
    w3 = v3 + gamma*y;
    
    R = sqrt(w1.^2 + w2.^2);
    I = R > beta;
    w1(I) = beta*w1(I)./R(I);
    w2(I) = beta*w2(I)./R(I);
    w3(w3 > 0) = 0;
    ATwold = ATw;
    ATw = div(w1,w2) + w3;
    
    y   = z - alphaX.*ATw;
    gstarv = 0;
    Psi    = -(D(:)'*y(:).^2)/2/alpha - gstarv + Psifixedpart;
    psivec(iter+1) = -Psi;
    ybar   = max(y,0);
    g1ybar = grad1(ybar); g2ybar = grad2(ybar); normgybar = sqrt(g1ybar.^2+g2ybar.^2);
    f1ybar = beta*sum(normgybar(:));
    hsigmaplus = (D(:)'*(ybar(:)-z(:)).^2)/2/alpha + f1ybar + hsigmafixedpart;
    hsigmavec(iter+1) = hsigmaplus;
    if verbose,
        fprintf('\n    %3g)  hsigmaplus %g psi %g',iter,hsigmaplus,Psi);
    end
    if hsigmaplus <= eta*Psi
        break;
    end

    %extrapolation step
    alphan = (iter -1)/(iter + a);
    v1 = w1 + alphan*(w1-w1old);
    v2 = w2 + alphan*(w2-w2old);
    v3 = w3 + alphan*(w3-w3old);
    w1old = w1; w2old = w2; w3old = w3;
    
    ATv    = ATw + alphan*(ATw - ATwold);%div(v1,v2)+ v3;
    y      = z - alphaX.*ATv;
    g1y    = grad1(y); g2y = grad2(y); 
    
end
G = hsigmaplus - Psi;
if nargout >= 11 
   psivec(iter +2:end) = [];
   hsigmavec(iter+2:end) = [];
   varargout{1}.psi = psivec; varargout{1}.hsigma = hsigmavec; 
end
