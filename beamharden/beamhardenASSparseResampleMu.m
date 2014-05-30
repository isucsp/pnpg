
function out = beamhardenASSparseResampleMu(Phi,Phit,Psi,Psit,y,xInit,opt)
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
%   Note: After the original algorithm converges, resample mu to get better 
%   reconstruction.
%
%   Reference:
%   Author: Renliang Gu (renliang@iastate.edu)
%   $Revision: 0.1 $ $Date: Fri 12 Jul 2013 11:49:11 AM CDT

maxItr=opt.maxItr;
opt.maxItr=floor(maxItr/opt.numCall);
out=[];

for i=1:opt.numCall
    tempOut=beamhardenASSparse(Phi,Phit,Psi,Psit,y,xInit,opt);
    out=[out; tempOut];

    %figure; semilogx(tempOut.mu, tempOut.Ie,'r*-');
    outIe=tempOut.Ie; outMu=log(tempOut.mu); xInit=tempOut.alpha;

    outIe=[outIe(:)'; zeros(1,length(outIe))]; outIe=outIe(:);
    outIe=[0; 0; outIe(1:end-1); 0; 0];
    temp=find(conv(outIe,ones(3,1),'same'));

    outMu=[outMu(:)'; zeros(1,length(outMu))]; outMu=outMu(1:end-1);
    outMu=conv(outMu,[0.5 1 0.5],'same');
    outMu=[outMu(1)*2-outMu(2); outMu(:); outMu(end)*2-outMu(end-1)];
    outMu=[outMu(1)*2-outMu(2); outMu(:); outMu(end)*2-outMu(end-1)];

    opt.mu=exp(outMu(temp)); opt.Ie=outIe(temp);
    %hold on; plot(opt.mu,opt.Ie,'+-');
    opt.E=length(opt.Ie);
    opt.sampleMode='assigned';
end

end
