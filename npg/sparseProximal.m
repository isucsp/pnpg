function proximalOut=sparseProximal(Psi, Psit, prj_C, method, opt)
% make an object to solve 0.5*||x-a||_2^2+u*||Psi^T(x)||_1+I_C(x)
% Psi, Psit and C are specified via constructor; see denoise for a and u

% default to have not convex set contraint C
if(~exist('prj_C','var') || isempty(prj_C)) prj_C=@(x)x; end
if(~exist('method','var') || isempty(method)) method='admm'; end
proxOp=@(p) min(max(p,-1),1);

% Note that in special cases, such as DWT, initial step can
% have better values to start from.
% In that case, set opt.initStep='fixed' with the approximated
% Lipschitz: opt.Lip=@(u)u^2; if Psi(Psit(x))==x.

proximalOut.val = @(x) norm(reshape(Psit(x),[],1),1);
proximalOut.exact=false;
switch(lower(method))
    case 'pnpg'
        proximalOut.prox=@denoisePNPG;
    case 'admm'
        proximalOut.prox=@denoiseADMM;
    otherwise
        proximalOut.prox=@denoiseADMM;
end

if(~exist('opt','var') || ~isfield(opt,'initStep')) opt.initStep='bb'; end
if(~isfield(opt,'maxLineSearch')) opt.maxLineSearch=5; end
if(~isfield(opt,'minItr')) opt.minItr=1; end

if(~isfield(opt,'stepIncre')) opt.stepIncre=0.5^0.2; end
if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
% By default disabled.  Remember to use a value around 5 for the Poisson model
% with poor initialization.
if(~isfield(opt,'preSteps')) opt.preSteps=0; end
if(~isfield(opt,'gamma')) gamma=2; else gamma=opt.gamma; end
if(~isfield(opt,'b')) b=0.25; else b=opt.b; end
if(~isfield(opt,'Lip')) opt.Lip=[]; end

if(~isfield(opt,'cumuTol')) opt.cumuTol=4; end
if(~isfield(opt,'incCumuTol')) opt.incCumuTol=true; end
if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=true; end
if(~isfield(opt,'backtracking')) opt.backtracking=true; end
if(~isfield(opt,'usePInit')) opt.usePInit=true; end
if(~isfield(opt,'dualGap')) opt.dualGap=false; end

if(opt.adaptiveStep) opt.backtracking=true; end

% Debug output information
% >=0: no print,
% >=1: only report results,
% >=2: detail output report, 
% >=4: plot real time cost and RMSE,
if(~isfield(opt,'debugLevel')) opt.debugLevel=0; end
% print iterations every opt.verbose lines.
if(~isfield(opt,'verbose')) opt.verbose=100; end

% Output options and debug information
% >=0: minimum output with only results,
% >=1: some cheap output,
% >=2: more detail output and expansive (impairs CPU time, only for debug)
if(~isfield(opt,'outLevel')) opt.outLevel=0; end

debug=Debug(opt.debugLevel);

function [x,itr,p,out]=denoisePNPG(a,u,thresh,maxItr,pInit)
    % For Psi and Psit that are from wavelet transform, set
    % opt.Lip=@(u)u^2;
    %
    if(isempty(opt.Lip) || ~isa(opt.Lip, 'function_handle'))
        Lip=inf;
    else
        Lip=opt.Lip(u);
    end

    if(nargout<4)
        opt.outLevel=0;
    end

    % print start information
    if(debug.level(2))
        fprintf('\n%s\n', repmat( '=', 1, 80 ) );
        str=sprintf('Projected Nestrov''s Proximal-Gradient (PNPG) Method');
        fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
        fprintf('%s\n', repmat('=',1,80));
        str=sprintf( ' %5s','Itr');
        str=sprintf([str ' %14s'],'Objective');
        str=sprintf([str ' %12s %4s'], '|dx|/|x|', 'lSrh');
        str=sprintf([str ' %12s'], '|d Obj/Obj|');
        if(opt.dualGap)
            str=sprintf([str ' %12s'], 'dualGap');
        end
        str=sprintf([str '\t u=%g'],u);
        fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
    end

    tStart=tic;

    if(~isempty(pInit) && iscell(pInit) && opt.usePInit)
        p=pInit; 
    else
        p=zeros(size(Psit(a)));
    end

    itr=0; convThresh=0; theta=1; preP=p;
    Psi_p=Psi(p); prePsi_p=Psi_p;
    y=a-u*Psi_p; % is real
    x=prj_C(y);   % is real
    cost=(sqrNorm(y)-sqrNorm(x-y))/2;
    goodStep=true;
    t=stepSizeInit(opt.initStep,Lip);

    if(opt.outLevel>=1) out.debug={}; end
    if(opt.adaptiveStep) cumu=0; end

    while(true)
        if(itr >= maxItr || convThresh>=1) break; end
        itr=itr+1;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  start of one PNPG step  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        numLineSearch=0; incStep=false; goodMM=true;
        if(opt.adaptiveStep)
            if(cumu>=opt.cumuTol)
                % adaptively increase the step size
                t=t*opt.stepIncre;
                cumu=0;
                incStep=true;
            end
        else
            if(itr==1)
                newTheta=1;
            else
                B=t/preT;
                newTheta=1/gamma+sqrt(b+B*theta^2);
                %newTheta=(1+sqrt(1+4*B*theta^2))/2;
            end
            pbar=p+(theta-1)/newTheta*(p-preP);
            ybar=a-u*(Psi_p  + (theta-1)/newTheta*(Psi_p-prePsi_p));
            xbar=prj_C(ybar);   % is real
            if(opt.backtracking)
                oldCost=(sqrNorm(ybar)-sqrNorm(xbar-ybar))/2;
            end
            grad=-u*Psit(xbar);
        end

        % start of line Search
        while(true)
            numLineSearch = numLineSearch+1;

            if(opt.adaptiveStep)
                if(itr==1)
                    newTheta=1;
                else
                    B=t/preT;
                    newTheta=1/gamma+sqrt(b+B*theta^2);
                    %newTheta=(1+sqrt(1+4*B*theta^2))/2;
                end
                pbar=p+(theta-1)/newTheta*(p-preP);
                ybar=a-u*(Psi_p  + (theta-1)/newTheta*(Psi_p-prePsi_p));
                xbar=prj_C(ybar);   % is real
                if(opt.backtracking)
                    oldCost=(sqrNorm(ybar)-sqrNorm(xbar-ybar))/2;
                end
                grad=-u*Psit(xbar);
            end

            newP=proxOp(pbar-grad/t);
            newPsi_p=Psi(newP);
            newY=a-u*newPsi_p; % is real
            newX=prj_C(newY);   % is real
            newCost=(sqrNorm(newY)-sqrNorm(newX-newY))/2;
            if(~opt.backtracking) break; end
            if((newCost-oldCost)<=...
                    innerProd(newP-pbar,grad)...
                    +t*sqrNorm(newP-pbar)/2)
                if(itr<=opt.preSteps && opt.adaptiveStep && goodStep)
                    cumu=opt.cumuTol;
                end
                break;
            else
                if(numLineSearch<=opt.maxLineSearch && t<Lip)
                    t=t/opt.stepShrnk; goodStep=false;
                    % Penalize if there is a step size increase just now
                    if(incStep)
                        incStep=false;
                        if(opt.incCumuTol)
                            opt.cumuTol=opt.cumuTol+4;
                        end
                    end
                else  % don't know what to do, mark on debug and break
                    goodMM=false;
                    debug.appendLog('_FalseMM');
                    break;
                end
            end
        end

        if((newCost-cost)>1e-14*norm([newCost,cost],inf))
            if(goodMM)
                if(restart()) itr=itr-1; continue; end
            end

            % give up and force it to converge
            debug.appendLog('_ForceConverge');
            preP=p; difX=0;
            prePsi_p=Psi_p;
            preCost=cost;
        else
            difX = relativeDif(x,newX);
            x=newX;
            preP = p; prePsi_p=Psi_p;
            p = newP; Psi_p=newPsi_p;
            theta = newTheta;
            preCost=cost;
            cost = newCost;
        end
        preT=t;

        if(opt.dualGap)
            % slightly larger than the exact gap, but saves time
            dualGap=norm(grad(:),1)-u*sum(reshape(x.*Psi_p,[],1));
        end

        if(opt.adaptiveStep)
            if(numLineSearch==1)
                cumu=cumu+1;
            else
                cumu=0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %  end of one PNPG step  %
        %%%%%%%%%%%%%%%%%%%%%%%%%%

        if(opt.dualGap)
            if(dualGap<=thresh && itr>=opt.minItr)
                convThresh=convThresh+1;
            end
        else
            if(difX<=thresh && itr>=opt.minItr)
                convThresh=convThresh+1;
            end
        end

        if(opt.outLevel<1 && opt.debugLevel<2)
            continue;
        end

        difCost=abs(cost-preCost)/max(1,abs(cost));
        if(opt.outLevel>=1)
            out.time(itr)=toc(tStart);
            out.cost(itr)=cost;
            out.difX(itr)=difX;
            if(opt.dualGap)
                out.dualGap(itr)=dualGap;
            end
            out.difCost(itr)=difCost;
            out.theta(itr)=theta;
            out.numLineSearch(itr) = numLineSearch;
            out.stepSize(itr) = 1/t;
            if(~isempty(debug.log()))
                out.debug{size(out.debug,1)+1,1}=itr;
                out.debug{size(out.debug,1),2}=debug.log();
                debug.clearLog();
            end;
        end

        if(debug.level(2))
            debug.print(2,sprintf(' %5d',itr));
            debug.print(2,sprintf(' %14.8g',cost));
            debug.print(2,sprintf(' %12g %4d',difX,numLineSearch));
            debug.print(2,sprintf(' %12g', difCost));
            if(opt.dualGap)
                debug.print(2,sprintf(' %12g', dualGap));
            end
            if(~debug.clear_print(2))
                if(mod(itr,opt.verbose)==0) debug.println(2); end
            end
        end
    end

    if(nargout>=4)
        out.opt = opt;
        out.date=datestr(now);
        out.gap=norm(grad(:),1)-u*sum(reshape(x.*Psi_p,[],1));
    end
    if(opt.outLevel>=2)
        out.grad=grad;
    end
    if(debug.level(1))
        fprintf('\nCPU Time: %g, objective=%g\n',toc(tStart),cost);
    end

    function res=restart()
        % if has monmentum term, restart
        res=(pNorm(xbar-x,0)~=0);
        if(~res) return; end

        theta=1;
        debug.appendLog('_Restart');
        debug.printWithoutDel(2,'\t restart');
    end

    function t_=stepSizeInit(select,Lip_)
        switch (lower(select))
            case 'bb'   % use BB method to guess the initial stepSize
                y_=a-u*Psi(p); % is real
                x_=prj_C(y_);   % is real
                grad_=-u*Psit(x_);
                grad_=grad_/pNorm(grad_,2);
                t_ = u^2*sqrNorm(Psi(grad_));
            case 'fixed'
                t_ = Lip_;
            otherwise
                error('unkown selection for initial step: %s', select);
        end
        if(isnan(t_) || t_<=0)
            error('PNPG is having a negative or NaN step size, do nothing and return!!');
        end
    end
end % function [x,itr,p,out]=denoisePNPG(a,u,thresh,maxItr,pInit)


function [x,itr,p]=denoiseADMM(a,u,relativeTol,maxItr,pInit)
    %
    % solve 0.5*||α-a||_2^2 + I(α≥0) + u*||Psit(α)||_1
    %
    % author: Renliang Gu (gurenliang@gmail.com)
    %
    if((~exist('relativeTol','var')) || isempty(relativeTol)) relativeTol=1e-6; end
    if((~exist('maxItr','var')) || isempty(maxItr)) maxItr=1e3; end
    if((~exist('pInit','var')) || isempty(pInit) || ~iscell(pInit))
        temp=size(Psit(a));
        pInit={zeros(temp), zeros(temp), 1};
    end
    % this makes sure the convergence criteria is nontrival
    relativeTol=min(1e-3,relativeTol); strlen=0;

    % scale the input to prevent numerical problem
    scale=pNorm(a,2);
    if(scale==0) x=zeros(size(a)); itr=0; p=pInit; return; end
    s=pInit{1}; nu=pInit{2}; rho=pInit{3};

    a=a/scale;  u=u/scale; s=s/scale; nu=nu/scale;

    itr=0; cnt=0; preS=s;
    while(true)
        itr=itr+1;
        cnt= cnt + 1;

        x = prj_C(scale*(a+rho*Psi(s+nu))/(1+rho))/scale;
        Psit_x=Psit(x);
        s = Utils.softThresh(Psit_x-nu,u/rho);
        nu=nu+s-Psit_x;

        difS=pNorm(s-preS,2); preS=s;
        residual = pNorm(s-Psit_x,2);

        if(debug.level(2))
            cost=0.5*sqrNorm(max(Psi(s),0)-a)+u*pNorm(Psit(max(Psi(s),0)),1);
            gap=rho*nu'*(s-Psit_x);

            str=sprintf('itr=%d, cost=%g pRes=%g dRes=%g gap=%g rho=%g       ',itr,...
                cost,residual,difS,gap,rho);
            if(strlen==0 || (mod(itr-1,100)==0 || (itr<=100 && mod(itr-1,10)==0) || itr-1<10))
                fprintf('\n%s',str);
            else
                fprintf([repmat('\b',1,strlen) '%s'],str);
            end
            strlen = length(str);
        end

        if(itr>maxItr) break; end
        if(difS<=relativeTol && residual<=relativeTol) break; end
        if(cnt>10) % prevent excessive back and forth adjusting
            if(difS>10*residual)
                rho = rho/2 ; nu=nu*2; cnt=0;
            elseif(difS<residual/10)
                rho = rho*2 ; nu=nu/2; cnt=0;
            end
        end
    end 
    x = prj_C(scale*(a+rho*Psi(s+nu))/(1+rho));
    p = {s,nu,rho};
    % end of the ADMM inside the NPG
    if(debug.level(2))
        fprintf('\n');
    end
end % function [x,itr,p]=denoiseADMM(a,u,relativeTol,maxItr,pInit)

end % function proximalOut=sparseProximal(Psi, Psit, prj_C, method, opt)
