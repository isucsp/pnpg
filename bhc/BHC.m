classdef BHC < handle
    methods(Static)
        function out = NPG2(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NPG'; opt.IeStep='NPG'; opt.errorType=0;
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = NPG_AS(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NPG'; opt.IeStep='ActiveSet'; opt.errorType=0;
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end

        % nonlinear conjugate gradient method (Polak-Ribière)
        function out = NCG_PR_AS(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NCG_PR'; opt.errorType=0;
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end

        % skip Ie step by assuming Ie is known
        function out = NPG(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NPG'; opt.skipIe=true; opt.errorType=0;
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end

%       opt.alphaStep='NPG'; opt.adaptiveStep=100; opt.cumuTol=10;
%       opt.alphaStep='NPG'; opt.CenterB=false; 
%       opt.alphaStep='NPG'; opt.CenterB=true; opt.correctCenterB=false;
%       opt.alphaStep='NPG'; opt.CenterB=true; opt.correctCenterB=true;

    end
end
