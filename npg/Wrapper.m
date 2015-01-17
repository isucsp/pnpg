classdef Wrapper < handle
    methods(Static)
        function out = NPG(Phi,Phit,Psi,Psit,y,xInit,opt)
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
            A = @(xx) Phi(Psi(xx)); At = @(yy) Psit(Phit(yy));
            B = @(xx) xx;
            opt.trueAlpha=Psit(opt.trueAlpha);
            out=lasso(A,At,B,B,y,Psit(xInit),opt);
            out.alpha=Psi(out.alpha);
        end
        function out = FPC(Phi,Phit,Psi,Psit,y,xInit,opt)
            if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
            if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
            A = @(xx) Phi(Psi(xx)); At = @(yy) Psit(Phit(yy));
            AO= @(xx,mode) AA(xx,mode,A,At);
            option.x0=Psit(xInit);
            if(isfield(opt,'trueAlpha'))
                option.xs=Psit(opt.trueAlpha);
            end
            option.xtol=opt.thresh;
            option.mxitr=opt.maxItr;
            option.scale=false;
            out = fpc_bb_mod(length(option.x0),AO,y,1/opt.u,[],option);
            out.alpha = Psi(out.x);
            out.difAlpha = out.step;
            if(isfield(opt,'trueAlpha'))
                out.RMSE=out.n2re.^2;
                %sqrNorm(out.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
            end
            out.opt = opt;
            fprintf('fpc cost=%g, RMSE=%g\n',out.cost(end),out.RMSE(end));
            function out = AA(xxx, mode , Phi, Phit)
                if(mode==1) out = Phi(xxx); else out = Phit(xxx); end
            end
        end
        function out = FPCas(Phi,Phit,Psi,Psit,y,xInit,opt)
            if(~isfield(opt,'maxItr')) opt.maxItr=2e3; end
            if(~isfield(opt,'thresh')) opt.thresh=1e-6; end
            if(isfield(opt,'trueAlpha')) option.trueAlpha=Psit(opt.trueAlpha); end
            A = @(xx) Phi(Psi(xx)); At = @(yy) Psit(Phit(yy));
            AO=A_operator(A,At);
            option.x0=Psit(xInit);
            option.minK=floor(0.01*length(option.x0(:)));
            option.mxitr=opt.maxItr;
            option.gtol = 0;
            option.gtol_scale_x = opt.thresh;
            [s, out] = FPC_AS_mod(length(option.x0),AO,y,opt.u,[],option);
            out.alpha = Psi(s);
            out.fVal=[0.5*sqrNorm(Phi(out.alpha)-y);...
                sqrNorm(out.alpha.*(out.alpha<0));...
                pNorm(Psit(out.alpha),1)];
            if(isfield(opt,'trueAlpha') && ~isfield(out,'RMSE'))
                out.RMSE=sqrNorm(out.alpha-opt.trueAlpha)/sqrNorm(opt.trueAlpha);
            end
            if(~isfield(out,'cost')) out.cost=out.f; end
            if(~isfield(out,'time')) out.time=out.cpu; end
            out.opt = opt;
            fprintf('fpcas cost=%g, RMSE=%g\n',out.f(end),out.RMSE(end));
        end
        function out = SPIRAL(Phi,Phit,Psi,Psit,y,xInit,opt)
            subtolerance=1e-5;
            if(~isfield(opt,'verbose')) opt.verbose=100; end
            if(~isfield(opt,'bb')) opt.bb=zeros(size(y)); end
            [out.alpha, out.p, out.cost, out.reconerror, out.time,out.difAlpha] = ...
                SPIRALTAP_mod(y,Phi,opt.u,'penalty','ONB',...
                'AT',Phit,'W',Psi,'WT',Psit,'noisetype',opt.noiseType,...
                'initialization',xInit,'maxiter',opt.maxItr,...
                'miniter',0,'stopcriterion',3,...
                'tolerance',opt.thresh,'truth',opt.trueAlpha,...
                'bb',opt.bb,...
                'subtolerance',subtolerance,'monotone',1,...
                'saveobjective',1,'savereconerror',1,'savecputime',1,...
                'reconerrortype',3,'savedifalpha',1,...
                'savesolutionpath',0,'verbose',opt.verbose);
            out.fVal=[0.5*sqrNorm(Phi(out.alpha)-y);...
                sqrNorm(out.alpha.*(out.alpha<0));...
                pNorm(Psit(out.alpha),1)];
            out.RMSE=out.reconerror;
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
            ppsi = @(xxx,uuu,thrsh) FISTA_ADMM_NNL1.adaptiveADMM(Psi,Psit,xxx,uuu,thrsh,100);
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
            % The following suggests that to get to the same RMSE, glmnet needs
            % far smaller convergence criteria
            options.thresh=opt.thresh;
            % options.maxit=opt.maxItr;

            switch lower(opt.noiseType)
                case lower('poissonLogLink')
                    L = @(aaa) Utils.poissonModelLogLink(aaa,@(xxx) Phi*xxx,@(xxx) Phi'*xxx,y);
                    [~,g]=L(xInit*0);
                    u_max=pNorm(Psi'*g,inf);
                    model='poisson';
                    options.intr = true; % need the interception
                    options.standardize=false;
                    A=-A;
                case 'gaussian'
                    L = @(aaa) Utils.linearModel(aaa,@(xxx) Phi*xxx,@(xxx) Phi'*xxx,y);
                    [~,g]=L(xInit*0);
                    u_max=pNorm(Psi'*g,inf);
                    u_max=pNorm(A'*y,inf);
                    options.intr = false;
                    options.standardize=false;
                    model='gaussian';
            end
            if(length(opt.u)>1)
                options.lambda=opt.u/length(y);
            end
            out=glmnet(A,y,model,options);

            out.time=toc;
            out.opt = opt;

            trueAlphaNorm=sqrNorm(opt.trueAlpha);
            for i=1:length(out.lambda)
                out.alpha(:,i) = Psi*out.beta(:,i);
                out.u(i)=out.lambda(i)*length(y);
                out.a(i)=log10(out.u(i)/u_max);
                if(isfield(opt,'trueAlpha'))
                    out.RMSE(i)=sqrNorm(out.alpha(:,i)-opt.trueAlpha)/trueAlphaNorm;
                end
                out.fVal(:,i)=[L(out.alpha(:,i));...
                    sqrNorm(out.alpha.*(out.alpha<0));...
                    pNorm(Psi'*(out.alpha),1)];
                out.cost(i)= out.fVal(1,i)+out.u(i)*out.fVal(3,i);
            end
            fprintf('glmnet cost=%g, RMSE=%g, time=%f\n',out.cost(end),min(out.RMSE()),out.time);
        end
    end
end

