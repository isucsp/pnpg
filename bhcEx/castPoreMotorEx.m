function castPoreMotorEx(op)
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
        clear('Oopt'); filename = [mfilename '.mat'];

        Oopt.beamharden=true; Oopt.spectBasis='b1'; Oopt.E=20;
        Oopt.estIe=true; Oopt.maxItr=4e3; Oopt.thresh=1e-6;

        prjFull = [60, 40, 72, 120, 180, 360];
        for i=length(prjFull):-1:1
            Oopt.prjFull = prjFull(i); Oopt.prjNum = Oopt.prjFull;

            [y,Phi,Phit,Psi,Psit,Oopt,FBP]=loadCastPoreMotor(Oopt);

            initSig = maskFunc(FBP(y),Oopt.mask~=0);

            j=1;
            fprintf('%s, i=%d, j=%d\n','Filtered Backprojection',i,j);
            fbp{i}.img=FBP(y);
            fbp{i}.alpha=fbp{i}.img(Oopt.mask~=0);

            if(i==6)
                j=3;
                u  =  10.^[-5  -5   -5   -5   -5   -5];
                opt=Oopt; opt.u=u(i)*10^(j-3); opt.proximal='tvl1';
                opt.alphaStep='NPG';
                opt.thresh=-1; opt.maxItr=1e4;
                npgTV_b1_long{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=Oopt; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.alphaStep='PG';
                opt.thresh=-1; opt.maxItr=1e4;
                pgTV_b1_long{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
                save(filename);
            end

            continue;

            % unknown ι(κ), NPG-AS
            for j=[3 4 2]
                fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
                %npg_b1{i,j}=BHC.NPG2(Phi,Phit,Psi,Psit,y,initSig,opt);
                u  =  10.^[-5  -5   -5   -5   -5   -5];
                opt=Oopt; opt.u=u(i)*10^(j-3); opt.proximal='tvl1';
                npgTV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

%               opt=Oopt; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.alphaStep='pg';
%               pgTV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=Oopt; opt.u=u(i)*10^(j-3); opt.proximal='wvltADMM';
                npgWV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

%               fpcas {i,j}=Wrapper.FPCas(Phi,Phit,Psi,Psit,y,initSig,opt);
                save(filename);
            end
        end

    case 'plot'
        load([mfilename '.mat']);
        prjFull = [60, 40, 72, 120, 180, 360];

        prjIdx=6; col=307; h=figure; forSave=[];

        img=showImgMask(      fbp{prjIdx     }.alpha,opt.mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'fbp_cpm.eps','psc2'); imwrite(img/maxImg,  'fbp_cpm.png');
        figure(h); plot(3*img(:,col),'b-'); hold on; forSave=[forSave, img(:,col)];

        aIdx=3; u  =  10.^[-5  -5  -5  -5  -5  -5];
        img=showImgMask(npgTV_b1{prjIdx,aIdx}.alpha,opt.mask); maxImg=min(4,max(img(:))); figure; showImg(img,0); saveas(gcf,'npgTV_cpm.eps','psc2'); imwrite(img/maxImg,'npgTV_cpm.png');
        fprintf('u for NPGTV is %e\n',npgTV_b1{prjIdx,aIdx}.opt.u);
        figure(h); plot(img(:,col),'g-.'); forSave=[forSave, img(:,col)];

        legend('FBP', 'NPG\_TV');
        save('profile_cpm.data','forSave','-ascii');

        clear('opt');
        a1=npgTV_b1{prjIdx,aIdx};
        [y,Phi,Phit,Psi,Psit,~,FBP]=loadCastPoreMotor(a1.opt);

        q=a1.kappa(2)/a1.kappa(1);
        polymodel=Spline(a1.opt.spectBasis,[a1.kappa(1)/q; a1.kappa(:); a1.kappa(end)*q]);
        polyIout = polymodel.polyIout;

        PhiAlpha=Phi(a1.alpha);
        PhiFbp=Phi(fbp{prjIdx}.alpha);

        s=linspace(min(PhiAlpha),max(PhiAlpha),100);
        idx=randi(length(PhiAlpha),1000,1);

        figure;
        plot(PhiAlpha(idx),y(idx),'.'); hold on;
        plot(PhiFbp(idx),y(idx),'g.');
        plot(s,-log(polyIout(s,a1.Ie)),'r-');
        legend('NPG-BFGS reconstruction', 'FBP reconstruction', 'fitted curve by NPG-BFGS');
        xlabel('\Phi\alpha');
        ylabel('I^{out}=-ln[ \int \iota(\kappa) exp( -\kappa\Phi\alpha ) d\kappa  ]');

        forSave=[PhiAlpha(idx),y(idx)];
        save('test1.data','forSave','-ascii');
        forSave=[PhiFbp(idx),y(idx)];
        save('test2.data','forSave','-ascii');
        forSave=[s(:), -log(polyIout(s,a1.Ie))];
        save('test3.data','forSave','-ascii');

        !for i in `seq 1 3`; do echo "" >> test$i.data; done
        !for i in `seq 1 3`; do echo "" >> test$i.data; done
        !cat test[1-3].data > linearization_castPoreMotor.data
        !rm test[1-3].data

        idx= find((PhiAlpha>3.8) & (PhiAlpha<4.2));
        figure; hist(exp(-y(idx)),100);

        keyboard

        paperDir = './';
        %system(['mv effectiveCenterB.data ' paperDir]);
end
end


