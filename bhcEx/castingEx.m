function castingEx(op)
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

        opt.beamharden=true; opt.spectBasis='b1'; opt.E=20;
        opt.estIe=true;

        prjFull = [60, 40, 72, 120, 180, 360];
        for i=6
            opt.prjFull = prjFull(i); opt.prjNum = opt.prjFull;

            [y,Phi,Phit,Psi,Psit,opt,FBP]=loadCasting(opt);
            opt.maxItr=4e3; opt.thresh=1e-6;

            initSig = maskFunc(FBP(y),opt.mask~=0);

            j=1;
            fprintf('%s, i=%d, j=%d\n','Filtered Backprojection',i,j);
            fbp{i}.img=FBP(y);
            fbp{i}.alpha=fbp{i}.img(opt.mask~=0);

            % unknown ι(κ), NPG-AS
            for j=[ 3]
                fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
                %npg_b1{i,j}=BHC.NPG2(Phi,Phit,Psi,Psit,y,initSig,opt);
                u  =  10.^[-5  -5   -5   -5   -5   -4.5];
                opt.u=u(i)*10^(j-3); opt.proximal='tvl1';
                % npgTV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
%               npgTV_b1{i,j}=BHC.NPG2(Phi,Phit,Psi,Psit,y,initSig,opt);

                u  =  10.^[-5  -5   -5   -5   -5   -4.5];
                opt.u=u(i)*10^(j-3); opt.proximal='wvltADMM';
%               npgWV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
                test=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

%               opt.alphaStep='NPGs';
%               npgsWV_dis{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

%               fpcas {i,j}=Wrapper.FPCas(Phi,Phit,Psi,Psit,y,initSig,opt);
                save(filename);
            end
        end

    case 'plot'
        load([mfilename '.mat']);
        prjFull = [60, 40, 72, 120, 180, 360];

        prjIdx=6; col=307; h=figure; forSave=[];

        img=showImgMask(      fbp{prjIdx     }.alpha,opt.mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'fbp_casting.eps','psc2'); imwrite(img/maxImg,  'fbp_casting.png');
        figure(h); plot(3*img(:,col),'b-'); hold on; forSave=[forSave, img(:,col)];

        aIdx=3; u  =  10.^[-5  -5  -5  -5  -5  -5];
        img=showImgMask(npgTV_b1{prjIdx,aIdx}.alpha,opt.mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'npgTV_casting.eps','psc2'); imwrite(img/maxImg,'npgTV_casting.png');
        fprintf('u for NPGTV is %e\n',10^(aIdx-3)*u(prjIdx));
        figure(h); plot(img(:,col),'g-.'); forSave=[forSave, img(:,col)];

        legend('FBP', 'NPG\_TV');
        save('profile_casting.data','forSave','-ascii');
        
        clear('opt');
        a1=npgTV_b1{prjIdx,aIdx};
        [y,Phi,Phit,Psi,Psit,~,FBP]=loadCasting(a1.opt);

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
        !cat test[1-3].data > linearization_casting.data
        !rm test[1-3].data

        idx= find((PhiAlpha>2.9) & (PhiAlpha<3.1));
        figure; hist(exp(-y(idx)),100);

        keyboard
        
        paperDir = './';
        %system(['mv effectiveCenterB.data ' paperDir]);
end
end


