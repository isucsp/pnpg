classdef NPG < handle
    methods(Static)
        function out = main(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='FISTA_ADMM_NNL1';
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = NPG_nads(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='FISTA_ADMM_NNL1'; opt.adaptiveStep=false;
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = PG(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='IST_ADMM_NNL1'; opt.adaptiveStep=true;
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = NPGc(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=true; opt.alphaStep='FISTA_ADMM_NNL1';
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = NPGc_nads(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=true; opt.alphaStep='FISTA_ADMM_NNL1'; opt.adaptiveStep=false;
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = PGc(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=true; opt.alphaStep='IST_ADMM_NNL1'; opt.adaptiveStep=true;
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = NPGsc(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=true; opt.alphaStep='FISTA_L1';
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = NPGs(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='FISTA_L1';
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = FISTA(Phi,Phit,Psi,Psit,y,xInit,opt)
            opt.continuation=false; opt.alphaStep='FISTA_L1'; opt.adaptiveStep=false;
            % to call FISTA, user need to specify the choice for the initial step size to be either 'fixed' or 'BB'
            % default is 'BB'
            % opt.initStep='fixed'; % 'BB'; %
            out=lasso(Phi,Phit,Psi,Psit,y,xInit,opt);
        end
        function out = FPCas(Phi,Phit,Psi,Psit,y,xInit,opt)
            if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
            if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
            A = @(xx) Phi(Psi(xx)); At = @(yy) Psit(Phit(yy));
            AO=A_operator(A,At);
            option.x0=Psit(xInit);
            option.mxitr=opt.maxItr;
            option.gtol = 1e-20; option.gtol_scale_x = opt.thresh;
            [s, out] = FPC_AS_mod(length(option.x0),AO,y,opt.u,[],option);
            out.alpha = Psi(s);
            out.fVal=[0.5*sqrNorm(Phi(out.alpha)-y);...
                sqrNorm(out.alpha.*(out.alpha<0));...
                pNorm(Psit(out.alpha),1)];
            if(isfield(opt,'trueAlpha'))
                out.RMSE=sqrNorm(out.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
            end
            out.opt = opt;
            fprintf('fpcas cost=%g, RMSE=%g\n',out.f,out.RMSE);
        end
        function out = SPIRAL(Phi,Phit,Psi,Psit,y,xInit,opt)
            subtolerance=1e-5;
            if(~isfield(opt,'verbose')) opt.verbose=100; end
            [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                SPIRALTAP_mod(y,Phi,opt.u,'penalty','ONB',...
                'AT',Phit,'W',Psi,'WT',Psit,'noisetype',opt.noiseType,...
                'initialization',xInit,'maxiter',opt.maxItr,...
                'miniter',0,'stopcriterion',3,...
                'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                'subtolerance',subtolerance,'monotone',1,...
                'saveobjective',1,'savereconerror',1,'savecputime',1,...
                'reconerrortype',3,'savedifalpha',1,...
                'savesolutionpath',0,'verbose',opt.verbose);
            out.fVal=[0.5*sqrNorm(Phi(out.alpha)-y);...
                sqrNorm(out.alpha.*(out.alpha<0));...
                pNorm(Psit(out.alpha),1)];
            out.opt=opt;
            fprintf('SPIRAL cost=%g, RMSE=%g\n',out.cost(end),out.reconerror(end));
        end

        function out = SpaRSA(Phi,Phit,Psi,Psit,y,xInit,opt)
            fprintf('SpaRSA start\n');
            if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
            ppsi = @(xxx,uuu,thrsh) Psi(Utils.softThresh(Psit(xxx),uuu));
            rrrr = @(xxx) pNorm(Psit(xxx),1);
            [x_SpaRSA,x_debias_SpaRSA,obj_SpaRSA,times_SpaRSA,debias_start_SpaRSA,out]=...
                SpaRSA_mod(y,Phi,opt.u,...
                'AT',Phit,...
                'Psi',ppsi,...
                'Phi',rrrr,...
                'Initialization',xInit,...
                'StopCriterion',5,...
                'ToleranceA',opt.thresh, ...
                'True_x',opt.trueAlpha,...
                'BB_variant',1,...
                'Safeguard',1,...
                'Monotone',0,...
                'Continuation',1,...
                'Verbose',opt.debugLevel>0,...
                'MaxiterA',opt.maxItr);
            out.alpha=x_SpaRSA; out.cost=obj_SpaRSA; out.time=times_SpaRSA;
            out.RMSE=out.mses/sqrNorm(opt.trueAlpha)*length(opt.trueAlpha);
            fprintf('SpaRSA cost=%g, RMSE=%g\n',out.cost(end),out.RMSE(end));
        end

        function out = SpaRSAp(Phi,Phit,Psi,Psit,y,xInit,opt)
            fprintf('SpaRSA nonnegative start\n');
            if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
            ppsi = @(xxx,uuu,thrsh) FISTA_ADMM_NNL1.innerADMM_v4(Psi,Psit,xxx,uuu,thrsh);
            rrrr = @(xxx) pNorm(Psit(xxx),1);
            xInit(xInit<0)=0;
            [x_SpaRSA,x_debias_SpaRSA,obj_SpaRSA,times_SpaRSA,debias_start_SpaRSA,out]=...
                SpaRSA_mod(y,Phi,opt.u,...
                'AT',Phit,...
                'Psi',ppsi,...
                'Phi',rrrr,...
                'Initialization',xInit,...
                'StopCriterion',5,...
                'ToleranceA',opt.thresh, ...
                'True_x',opt.trueAlpha,...
                'BB_variant',1,...
                'Safeguard',1,...
                'Monotone',0,...
                'Continuation',1,...
                'Verbose',opt.debugLevel>0,...
                'MaxiterA',opt.maxItr);
            out.alpha=x_SpaRSA; out.cost=obj_SpaRSA; out.time=times_SpaRSA;
            out.RMSE=out.mses/sqrNorm(opt.trueAlpha)*length(opt.trueAlpha);
            fprintf('SpaRSA nonnegative cost=%g, RMSE=%g\n',out.cost(end),out.RMSE(end));
        end
        function out = glmnet(Phi,Psi,y,xInit,opt)
            if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
            if(~isfield(opt,'thresh')) opt.thresh=1e-6; end

            tic;
            A = Phi*Psi;
            options=glmnetSet;
            options.thresh=opt.thresh;
            options.maxit=opt.maxItr;
            if(length(opt.u)>1)
                options.nlambda=length(opt.u);
                options.lambda=opt.u/length(y);
            end

            switch lower(opt.noiseType)
                case lower('poissonLogLink')
                    out=glmnet(A,y,'poisson',options);
                case 'gaussian'
                    options.intr = false;
                    options.standardize=false;
                    out=glmnet(A,y,'gaussian',options);
            end

            out.time=toc;
            out.opt = opt;

            trueAlphaNorm=sqrNorm(opt.trueAlpha);
            for i=1:length(out.lambda)
                out.alpha(:,i) = Psi*out.beta(:,i);
                out.u(i)=out.lambda(i)*length(y);
                if(isfield(opt,'trueAlpha'))
                    out.RMSE(i)=sqrNorm(out.alpha(:,i)-opt.trueAlpha)/trueAlphaNorm;
                end
                switch lower(opt.noiseType)
                    case 'gaussian'
                        out.fVal(:,i)=[0.5*sqrNorm(Phi*out.alpha(:,i)-y);...
                            sqrNorm(out.alpha(:,i).*(out.alpha(:,i)<0));...
                            pNorm(Psi'*out.alpha(:,i)),1)];
                        out.cost(i)= out.fVal(1,i)+out.u(i)*out.fVal(3,i);
                    case lower('poissonLogLink')
                        PhiAlpha=Phi*out.alpha(:,i); weight = exp(-PhiAlpha);
                        sumy=sum(y); sumWeight=sum(weight);
                        f=sumy*log(sumWeight)+innerProd(y,PhiAlpha);
                        out.fVal=[f;...
                            sqrNorm(out.alpha.*(out.alpha<0));...
                            pNorm(Psit(out.alpha),1)];
                        out.cost(i)= out.fVal(1,i)+out.u(i)*out.fVal(3,i);
                end
            end
            fprintf('glmnet cost=%g, RMSE=%g, time=%f\n',out.f(end),out.RMSE(end),out.time);
        end
    end
end

