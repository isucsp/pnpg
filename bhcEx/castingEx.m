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
        clear('OPT'); filename = [mfilename '.mat'];

        OPT.beamharden=true; OPT.spectBasis='b1'; OPT.E=20;
        OPT.maxItr=4e3; OPT.thresh=1e-6;
        OPT.estIe=true;

        prjFull = [60, 40, 72, 120, 180, 360];
        u  =  10.^[-5  -5   -5   -5   -5   -4.5];
        for i=6:-1:4
            OPT.prjFull = prjFull(i); OPT.prjNum = OPT.prjFull;
            [y,Phi,Phit,Psi,Psit,OPT,FBP]=loadCasting(OPT);
            initSig = maskFunc(FBP(y),OPT.mask~=0);

            j=1;
            fprintf('%s, i=%d, j=%d\n','Filtered Backprojection',i,j);
            fbp{i}.img=FBP(y);
            fbp{i}.alpha=fbp{i}.img(OPT.mask~=0);

            keyboard
            % unknown ι(κ), NPG-LBFGSB without sparsity constraints
            if(i==6)

                opt=OPT; opt.u=0; j=1; opt.alphaStep='NPG'; opt.proximal='nonneg';
                opt.trueAlpha=npgTV_b1{6,4}.alpha;
                npgTV_b1_u0{i,j}=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
                save(filename);
                continue;

                opt=OPT; opt.u=0; j=1; opt.alphaStep='NPG'; opt.proximal='nonneg';
                opt.maxItr=500; % this one works the best, see npgTV_b1_u0
                npgTV_b1_u0_i3=beamhardenSpline(Phi,Phit,Psi,Psit,y,initSig,opt);
                opt=OPT;
            end
            continue;


            for j=[5:-1:1]
                fprintf('%s, i=%d, j=%d\n','NPG-AS',i,j);
                %npg_b1{i,j}=BHC.NPG2(Phi,Phit,Psi,Psit,y,initSig,opt);

                if(i>=6 && j==3)
                    opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.alphaStep='NPG';
                    opt.thresh=0e-16; opt.maxItr=1e4;

                    npgTV_b1_long{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
                    save(filename);

                    opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.alphaStep='PG';
                    opt.thresh=0e-16; opt.maxItr=1e4;
                    pgTV_b1_long{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
                    save(filename);

%                   opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.restart=false;
%                   opt.thresh=1e-16; opt.maxItr=1e4;
%                   npgTV_b1_norestart{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
%                   save(filename);
                end

                continue;

                opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1';
                if(j<5)
                    opt.Ie=npgTV_b1_continuation{i,j+1}.Ie;
                    npgTV_b1_continuation{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,npgTV_b1_continuation{i,j+1}.alpha,opt);
                else
                    npgTV_b1_continuation{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
                end

                opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1';
%               npgTV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.saveAnimate=true;
                npgTV_b1_anim=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3+0.25); opt.proximal='tvl1';
                npgTV_b1_25{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3+0.5); opt.proximal='tvl1';
                npgTV_b1_50{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3+0.75); opt.proximal='tvl1';
                npgTV_b1_75{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='wvltADMM';
                npgWV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3+0.75); opt.proximal='wvltADMM';
                npgWV_b1_75{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

                opt=OPT; opt.u=u(i)*10^(j-3); opt.proximal='tvl1'; opt.alphaStep='PG';
                pgTV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);


%               opt.alphaStep='NPGs';
%               npgsWV_dis{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);

%               fpcas {i,j}=Wrapper.FPCas(Phi,Phit,Psi,Psit,y,initSig,opt);
                save(filename);
            end
        end

    case 'plot'
        load([mfilename '.mat']);
        prjFull = [60, 40, 72, 120, 180, 360];

        prjIdx=4; aIdx=4; row1=337; row2=531;
        h1=figure; forSave=[]; h2=figure;
        img=showImgMask(      fbp{prjIdx     }.alpha,opt.mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'fbp_casting_120.eps','psc2'); imwrite(img/maxImg,  'fbp_casting_120.png');
        figure(h1); plot(img(row1,:)/maxImg,'b-'); hold on; forSave=[forSave, reshape(img(row1,:),[],1)];
        figure(h2); plot(img(row2,:)/maxImg,'b-'); hold on; forSave=[forSave, reshape(img(row2,:),[],1)];

        img=showImgMask(npgTV_b1_25{prjIdx,aIdx}.alpha,opt.mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,'npgTV_casting_120.eps','psc2'); imwrite(img/maxImg,'npgTV_casting_120.png');
        fprintf('u for NPGTV is %e\n',npgTV_b1_25{prjIdx,aIdx}.opt.u);
        figure(h1); plot(img(row1,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row1,:),[],1)];
        figure(h2); plot(img(row2,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row2,:),[],1)];
        legend('FBP', 'NPG\_TV');
        save('profile_casting_120.data','forSave','-ascii');


        prjIdx=6; h1=figure; forSave=[]; h2=figure;
        img=showImgMask(      fbp{prjIdx     }.alpha,opt.mask); maxImg=max(img(:)); figure; showImg(img,0); saveas(gcf,  'fbp_casting.eps','psc2'); imwrite(img/maxImg,  'fbp_casting.png');
        figure(h1); plot(img(row1,:)/maxImg,'b-'); hold on; forSave=[forSave, reshape(img(row1,:),[],1)];
        figure(h2); plot(img(row2,:)/maxImg,'b-'); hold on; forSave=[forSave, reshape(img(row2,:),[],1)];

        aIdx=3; u  =  10.^[-5  -5  -5  -5  -5  -5];
        img=showImgMask(npgTV_b1{prjIdx,aIdx}.alpha,opt.mask); maxImg=max(img(:));
        fprintf('u for NPGTV is %e\n',10^(aIdx-3)*u(prjIdx)); temp=img(:)/maxImg;
        figure; showImg(img,0); saveas(gcf,'npgTV_casting.eps','psc2'); imwrite(img/maxImg,'npgTV_casting.png');
        figure(h1); plot(img(row1,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row1,:),[],1)];
        figure(h2); plot(img(row2,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row2,:),[],1)];

        img =showImgMask(npgTV_b1{prjIdx,aIdx+1}.alpha,opt.mask); maxImg=(img(:)'*img(:))/(img(:)'*temp);
        figure; showImg(img,0,maxImg);  imwrite(img/maxImg,     'npgTV_casting_u-4.png');
        figure(h1); plot(img(row1,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row1,:),[],1)];
        figure(h2); plot(img(row2,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row2,:),[],1)];

        img =showImgMask(npgTV_b1_u0{prjIdx}.alpha,opt.mask); maxImg=(img(:)'*img(:))/(img(:)'*temp);
        figure; showImg(img,0,maxImg);  imwrite(img/maxImg,     'npgTV_casting_u0.png');
        figure(h1); plot(img(row1,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row1,:),[],1)];
        figure(h2); plot(img(row2,:)/maxImg,'g-.'); forSave=[forSave, reshape(img(row2,:),[],1)];

        figure(h1); legend('FBP', 'NPG\_TV','NPG\_TV\_u=1e-4','NPG\_TV\_u=0');
        figure(h2); legend('FBP', 'NPG\_TV','NPG\_TV\_u=1e-4','NPG\_TV\_u=0');
        save('profile_casting.data','forSave','-ascii');

        keyboard
        
        clear('opt');
        a1=npgTV_b1{prjIdx,aIdx};
        a2=pgTV_b1{prjIdx,aIdx};

        forSave=[]; t=0;
        t=t+1; temp=a1.stepSize(:); forSave(1:length(temp),t)=temp;
        t=t+1; temp=a1.cost(:);     forSave(1:length(temp),t)=temp;
        t=t+1; temp=a1.time(:);     forSave(1:length(temp),t)=temp;
        t=t+1; temp=a2.stepSize(:); forSave(1:length(temp),t)=temp;
        t=t+1; temp=a2.cost(:);     forSave(1:length(temp),t)=temp;
        t=t+1; temp=a2.time(:);     forSave(1:length(temp),t)=temp;
        save('castingTrace.data','forSave','-ascii');

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

        Imea = exp(-y(:));
        Iout = polyIout(PhiAlpha,a1.Ie);

        residual = (Imea-Iout)./sqrt(Iout);
        [npgHist,npgRes] = hist(residual,500);

        Iout = exp(-PhiFbp);
        residual = (Imea-Iout)./sqrt(max(Iout,eps));
        [fbpHist,fbpRes] = hist(residual,500);

        forSave=[npgRes,npgHist,fbpRes,fbpHist];
        save('castingRes.data','forSave','-ascii');

        figure; plot(npgRes,npgHist); hold on; plot(fbpRes,fbpHist,'r');
        legend('npg','fbp');

        keyboard
        
        paperDir = './';
        %system(['mv effectiveCenterB.data ' paperDir]);
end
end

function [] = selectRegPar()
    npgTV_b1{i,j}=BHC.main(Phi,Phit,Psi,Psit,y,initSig,opt);
end
