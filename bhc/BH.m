classdef BH < handle
    methods(Static)
        function out = NPG_AS(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NPG';
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end

        % nonlinear conjugate gradient method (Polak-Ribière)
        function out = NCG_PR_AS(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NCG_PR';
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end

        % skip Ie step by assuming Ie is known
        function out = NPG(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='NPG'; opt.skipIe=true;
            out=beamhardenSpline(Phi,Phit,Psi,Psit,y,xInit,opt);
        end

        opt.alphaStep='NPG'; opt.adaptiveStep=100; opt.cumuTol=10;
        opt.alphaStep='NPG'; opt.CenterB=false; 
        opt.alphaStep='NPG'; opt.CenterB=true; opt.correctCenterB=false;
        opt.alphaStep='NPG'; opt.CenterB=true; opt.correctCenterB=true;

    end
end
