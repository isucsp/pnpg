function slPoiLogEx(op)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Reconstruction of Nonnegative Sparse Signals Using Accelerated
%                      Proximal-Gradient Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Author: Renliang Gu (renliang@iastate.edu)
%   v_0.2:      Changed to class oriented for easy configuration
%
% Skyline log link Poisson example
% vary the number of measurements and inital intensity constatn I_0

if(~exist('op','var')) op='plot'; end

switch lower(op)
    case {lower('I016'), lower('I012'), lower('I009')}
        if(strcmpi(op,'I016'))
            filename = [mfilename '_007.mat'];
            if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end;
            clear('opt'); opt.B=2^16; aa=-3:-0.1:-6;
        elseif(strcmpi(op,'I012'))
            filename = [mfilename '_027.mat'];
            if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end;
            clear('opt'); opt.B=2^12; aa=-1:-0.2:-5;
        elseif(strcmpi(op,'I009'))
            filename = [mfilename '_037.mat'];
            if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end;
            clear('opt'); opt.B=2^9; aa=-0.2:-0.1:-4;
        end

        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e5; opt.thresh=1e-6;
        m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200
        a=-4.4:-0.2:-4.8;
        for k=1:1
            for i=1:length(m)
                opt.m=m(i); opt.noiseType='poissonLogLink'; opt.matrixType='conv';
                [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
                initSig=-Phit(invEAAt*log(max(y,1)/max(y)))*0;
                fprintf('i=%d, k=%d, min=%d, max=%d\n',i,k,min(y), max(y));

                %% fit by approximated Gaussian model for known I0
                % u_max=pNorm(Psit(Phit(y.*log(opt.I0./max(y,1)))),inf);
                %% fit by Poisson model with known I0
                % u_max=pNorm(Psit(Phit(y-opt.I0)),inf);
                %% fit by the unkown I0 log link Poisson model
                u_max=pNorm(Psit(Phit(y-mean(y))),inf);

                %               opt.noiseType='poissonLogLink';
                %               temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
                %               xx=opt.thresh; opt.thresh=1e-12;
                %               gnet    {i,k}=Wrapper.glmnet(Utils.getMat(Phi,length(initSig)),Utils.getMat(Psi,length(Psit(initSig))),y,initSig,opt);
                %               opt.thresh=xx; clear('xx');
                %               npgFull  {i,k}=Wrapper.NPG  (Phi,Phit,Psi,Psit,y,initSig,opt);
                %               npgsFull {i,k}=Wrapper.NPGs (Phi,Phit,Psi,Psit,y,initSig,opt);
                %               opt=temp;

                opt.noiseType='poissonLogLink0';
                temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
                %               npgFull_knownI0 {i,k}=Wrapper.NPG (Phi,Phit,Psi,Psit,y,initSig,opt);
                %               npgsFull_knownI0{i,k}=Wrapper.NPGs(Phi,Phit,Psi,Psit,y,initSig,opt);
                xx=opt.thresh; opt.thresh=1e-12;
                gnet0   {i,k}=Wrapper.glmnet(Utils.getMat(Phi,length(initSig)),Utils.getMat(Psi,length(Psit(initSig))),y,initSig,opt);
                opt.thresh=xx; clear('xx');
                opt=temp;

                %               opt.noiseType='gaussian';
                %               %% fit by approximated Gaussian model for known I0
                %               yy=log(opt.I0./max(y,1)); yy=yy.*sqrt(y);
                %               wPhi=@(xxx) sqrt(y).*Phi(xxx); wPhit=@(xxx) Phit(sqrt(y).*xxx);

                %               temp=opt; opt.fullcont=true; opt.u=10.^aa*u_max;
                %               npglwFull {i,k}=Wrapper.NPG(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
                %               npgslwFull{i,k}=Wrapper.NPGs(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
                %               opt=temp;

                save(filename);
            end
        end

    case lower('ind')

        filename = [mfilename '_017.mat'];
        if(~exist(filename,'file')) save(filename,'filename'); else load(filename); end
        clear('opt');
        RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',0));
        opt.maxItr=1e5; opt.thresh=1e-6;
        m=[ 200, 300, 400, 500, 600, 700, 800, 900, 1024]; % should go from 200
        a=-4.4:-0.2:-4.8;
        for k=1:5
            for i=1:length(m)
                opt.m=m(i); opt.noiseType='poissonLogLink'; opt.matrixType='conv';
                [y,Phi,Phit,Psi,Psit,opt,~,invEAAt]=loadLinear(opt);
                initSig=-Phit(invEAAt*log(max(y,1)/max(y)))*0;
                fprintf('i=%d, k=%d, min=%d, max=%d\n',i,k,min(y), max(y));

                %% fit by approximated Gaussian model for known I0
                % u_max=pNorm(Psit(Phit(y.*log(opt.I0./max(y,1)))),inf);
                yy=log(opt.I0./max(y,1)); yy=yy.*sqrt(y);
                wPhi=@(xxx) sqrt(y).*Phi(xxx); wPhit=@(xxx) Phit(sqrt(y).*xxx);
                %% fit by Poisson model with known I0
                % u_max=pNorm(Psit(Phit(y-opt.I0)),inf);
                %% fit by the unkown I0 log link Poisson model
                u_max=pNorm(Psit(Phit(y-mean(y))),inf);

                for j=1:length(a)
                    opt.u=10^a(j)*u_max;

                    opt.noiseType='poissonLogLink';
                    npg {i,j,k}=Wrapper.NPGc   (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npgs{i,j,k}=Wrapper.NPGsc  (Phi,Phit,Psi,Psit,y,initSig,opt);

                    opt.noiseType='poissonLogLink0';
                    npg0 {i,j,k}=Wrapper.NPGc (Phi,Phit,Psi,Psit,y,initSig,opt);
                    npgs0{i,j,k}=Wrapper.NPGsc(Phi,Phit,Psi,Psit,y,initSig,opt);

                    opt.noiseType='gaussian';
                    npglw {i,j,k}=Wrapper.NPGc(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
                    npgslw{i,j,k}=Wrapper.NPGsc(wPhi,wPhit,Psi,Psit,yy,initSig,opt);
                end
                save(filename);
            end
        end

    case 'plot'

    end
end

