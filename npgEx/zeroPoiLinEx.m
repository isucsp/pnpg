function zeroPoiLinEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration
%
%                    Zero signal Poisson linear Model Example
%    vary the number of measurements and inital intensity constatn I_0

% zero signal Poisson linear model example, no background noise
% vary the number of measurements, with continuation

opt.maxItr=1e4; opt.thresh=1e-6; opt.debugLevel=1;
opt.noiseType='poisson'; opt.matrixType='nonneg'; opt.snr=1e7;
m=[ 200, 400, 600, 800];

for i=1:length(m)
    opt.m=m(i);
    [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
    opt.trueAlpha=zeros(size(opt.trueAlpha)); y=y*0;
    fprintf('min=%d, max=%d, sum(y)=%d\n',min(y), max(y),sum(y));

    initSig=Phit(invEAAt*y)*0+1;
    initSig=max(randn(size(initSig)),0);
    opt.u=1e-33;
    opt.u=pNorm(Psit(Phit(ones(size(y)))),inf)*2;

    Tnpgs  =Wrapper.NPGs  (Phi,Phit,Psi,Psit,y,initSig,opt); norm(Tnpgs.alpha)
    Tnpg   =Wrapper.NPG   (Phi,Phit,Psi,Psit,y,initSig,opt); norm(Tnpg.alpha)
    %Tspiral=Wrapper.SPIRAL(Phi,Phit,Psi,Psit,y,initSig,opt); norm(Tspiral.alpha)

end

end
