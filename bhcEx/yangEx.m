function yangEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polychromatic Sparse Image Reconstruction and Mass Attenuation Spectrum 
%            Estimation via B-Spline Basis Function Expansion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:  Changed to class oriented for easy configuration

% dis discretized, polychromatic model
% dis, compare the NPG and NCG_PR for both continuation and 
% non-continuation
%
% compare the effect forcing center model of spectrum 
%
% A few thing to tune: continuation, centerb, ActiveSet VS FISTA_Simplex
% maxIeSteps

if(~exist('op','var')) op='run'; end

switch lower(op)
    case 'run'
        filename = [mfilename '.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt'); filename = [mfilename '.mat'];

        opt.beamharden=true; opt.spectBasis='b1'; opt.E=30;
        opt.estIe=true;

        prjFull = [60, 80, 100, 120, 180, 360];
        for i=[1,length(prjFull)]
            opt.prjFull = prjFull(i); opt.prjNum = opt.prjFull;

            [y,Phi,Phit,Psi,Psit,opt,FBP]=loadYang(opt);
            opt.maxItr=4e3; opt.thresh=1e-6;

            initSig = maskFunc(FBP(y),opt.mask~=0);

            j=1;
            fprintf('%s, i=%d, j=%d\n','Filtered Backprojection',i,j);
            fbp{i}.img=FBP(y);
            fbp{i}.alpha=fbp{i}.img(opt.mask~=0);
            fbp{i}.RMSE=1-(innerProd(fbp{i}.alpha,opt.trueAlpha)^2)/sqrNorm(opt.trueAlpha)/sqrNorm(fbp{i}.alpha);
            fprintf('fbp RMSE=%g\n',fbp{i}.RMSE);
             
            % known ι(κ), NPG
            for j=1:5
                fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
                u  =  10.^[-5  -5   -5   -5   -5   -5];
                opt.u=10^(j-3)*u(i);
                opt.alphaStep='NPG'; opt.proximal='tviso'; opt.skipIe=true;
                npgTValpha_b1{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.skipIe=false;
%               opt.alphaStep='NPG'; opt.proximal='wvltADMM'; opt.skipIe=true;
%               npgWValpha_b1{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
%               opt.skipIe=false;
            end

            save(filename);
            continue

            % known ι(κ), linearization
            kappa=logspace(-floor(opt.E/2)/(opt.E-1)*3,...
                floor(opt.E/2-0.5)/(opt.E-1)*3,opt.E);
            q=kappa(2)/kappa(1);
            polymodel=Spline(opt.spectBasis,[kappa(1)/q; kappa(:); kappa(end)*q]);
            polyIout = polymodel.polyIout; clear('q');
            [opt.upkappa,opt.upiota]=getUpiota(opt.epsilon,opt.kappa,opt.iota);
            trueIe=interp1(log(opt.upkappa), opt.upiota ,log(kappa(:)),'spline');
            trueIe=max(0,trueIe);
            s=linspace(min(y(:))/10,max(y(:))*10,10000);
            yy=interp1(-log(polyIout(s,trueIe)),s,y,'spline');

            fprintf('%s, i=%d,\n','Linearized Filtered Backprojection',i);
            linFbp{i}.img=FBP(yy);
            linFbp{i}.alpha=linFbp{i}.img(opt.mask~=0);
            linFbp{i}.RMSE=1-(innerProd(linFbp{i}.alpha,opt.trueAlpha)^2)/sqrNorm(opt.trueAlpha)/sqrNorm(linFbp{i}.alpha);
            fprintf('linFbp RMSE=%g\n',linFbp{i}.RMSE);
            initSig = maskFunc(FBP(yy),opt.mask~=0);

            for j=4
                u = 10.^[-5  -5   -5   -5   -5   -5];
                opt.u=10^(j-3)*u(i)*max(abs(Psit(Phit(yy)))); opt.proximal='tviso';
                linNpgTV{i,j}=Wrapper.NPGc(Phi,Phit,Psi,Psit,yy,initSig,opt);
            end
             
            % unknown ι(κ), NPG-LBFGSB
            for j=[4]
                fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
                u  =  10.^[-5  -5   -5   -5   -5   -5];
                opt.u=10^(j-3)*u(i);
                npg_b1{i,j}=BHC.NPG2(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.alphaStep='NPG'; opt.proximal='tvl1';
                npgTV_b1{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.alphaStep='NPG'; opt.proximal='wvltADMM';
                npgWV_b1{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.alphaStep='NPGs'; opt.proximal='tvl1';
                npgsTV_b1{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.alphaStep='NPGs'; opt.proximal='wvltADMM';
                npgsWV_b1{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);

                save(filename);
            end

            % linear sparse model
            for j=1:5
                u  =  10.^[-5  -5   -5   -5   -5   -5];
                opt.u=10^(j-3)*u(i)*max(abs(Psit(Phit(yy))));
                opt.proximal='tvl1';
                npgTV{i,j}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt.proximal='wvltADMM';
                npgWV{i,j}=Wrapper.NPG(Phi,Phit,Psi,Psit,y,initSig,opt);
            end
        end

    case 'plot'
        load([mfilename '.mat']);


        
        paperDir = './';
        %system(['mv effectiveCenterB.data ' paperDir]);
end
end


