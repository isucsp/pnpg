function proximalOut=tvProximal(tv, prj_C)
    % make an object to solve 0.5*||x-a||_2^2+u*TV(x)+I_C(x)
    % TV and C are specified via constructor see denoise for a and u

    if(~exist('tv','var') || isempty(tv))
        tv='iso';
    end
    if(~exist('prj_C','var') || isempty(prj_C))
        prj_C=@(x)x;
    end
    switch(lower(tv))
        case 'iso'
            proxOp=@isoPrj;
        case 'l1'
            proxOp=@l1Prj;
    end

    proximalOut.penalty = @(z) tlv(z,tv);
    proximalOut.iterative=true;
    if(true)
        proximalOut.op=@denoisePNPG;
    else
        proximalOut.op=@denoise;
        prnt = 0;
    end

    if(~exist('opt','var') || ~isfield(opt,'adaptiveStep')) opt.adaptiveStep=true; end
    if(~isfield(opt,'debugLevel')) opt.debugLevel=0; end
    if(~isfield(opt,'outLevel')) opt.outLevel=0; end
    if(~isfield(opt,'initStep')) opt.initStep='fixed'; end
    if(~isfield(opt,'preSteps')) opt.preSteps=0; end
    if(~isfield(opt,'maxLineSearch')) opt.maxLineSearch=5; end
    if(~isfield(opt,'minItr')) opt.minItr=1; end

    function [x,itr,pOut,out]=denoisePNPG(a,u,thresh,maxItr,pInit)
        function [f,gp,gq,h] = dualFunc(p,q)
            % minimize 0.5*||a-u*real(Psi(p))||_2^2-0.5*||X-a+u*real(Psi(p))||_2^2
            % subject to ||p||_infty <= 1
            % where X=prj_C( a-u*real(Psi(p)) ), p and Psi may be complex.
            Qp=a-u*(Psi_v(p)+Psi_h(q)); % is real
            x=prj_C(Qp);   % is real
            f=(sqrNorm(Qp)-sqrNorm(x-Qp))/2;
            if(nargout>1)
                gp=-u*Psi_vt(x);
                gq=-u*Psi_ht(x);
                if(nargout>3)
                    h=[];
                end
            end
        end

        opt.Lip=8*u^2;

        if(isempty(pInit) || ~iscell(pInit))
            [I,J]=size(a);
            pInit={zeros(I-1,J); zeros(I,J-1)};
        end

        % default to not use any constraints.

        if(~exist('opt','var') || ~isfield(opt,'stepIncre')) opt.stepIncre=0.9; end
        if(~isfield(opt,'stepShrnk')) opt.stepShrnk=0.5; end
        % By default disabled.  Remember to use a value around 5 for the Poisson model
        % with poor initialization.
        if(~isfield(opt,'preSteps')) opt.preSteps=0; end
        if(~isfield(opt,'initStep')) opt.initStep='hessian'; end
        % Threshold for relative difference between two consecutive x
        if(~isfield(opt,'minItr')) opt.minItr=10; end
        if(~isfield(opt,'maxLineSearch')) opt.maxLineSearch=20; end
        if(~isfield(opt,'gamma')) gamma=2; else gamma=opt.gamma; end
        if(~isfield(opt,'b')) b=0.25; else b=opt.b; end

        if(~isfield(opt,'cumuTol')) opt.cumuTol=4; end
        if(~isfield(opt,'incCumuTol')) opt.incCumuTol=true; end
        if(~isfield(opt,'adaptiveStep')) opt.adaptiveStep=true; end

        % Debug output information
        % >=0: no print,
        % >=1: only report results,
        % >=2: detail output report, 
        % >=4: plot real time cost and RMSE,
        if(~isfield(opt,'debugLevel')) opt.debugLevel=1; end
        % print iterations every opt.verbose lines.
        if(~isfield(opt,'verbose')) opt.verbose=100; end

        % Output options and debug information
        % >=0: minimum output with only results,
        % >=1: some cheap output,
        % >=2: more detail output and expansive (impairs CPU time, only for debug)
        if(~isfield(opt,'outLevel')) opt.outLevel=0; end

        debug=Debug(opt.debugLevel);

        % print start information
        if(debug.level(2))
            fprintf('\n%s\n', repmat( '=', 1, 80 ) );
            str=sprintf('Projected Nestrov''s Proximal-Gradient (PNPG) Method');
            fprintf('%s%s\n',repmat(' ',1,floor(40-length(str)/2)),str);
            fprintf('%s\n', repmat('=',1,80));
            str=sprintf( ' %5s','Itr');
            str=sprintf([str ' %14s'],'Objective');
            str=sprintf([str ' %12s %4s'], '|dx|/|x|', 'αSrh');
            str=sprintf([str ' %12s'], '|d Obj/Obj|');
            str=sprintf([str '\t u=%g'],u);
            fprintf('%s\n%s\n',str,repmat( '-', 1, 80 ) );
        end

        tStart=tic;

        p=pInit{1}; q=pInit{2};
        itr=0; convThresh=0; theta=1; preP=p; preQ=q;
        cost=dualFunc(p,q);
        goodStep=true;
        t=stepSizeInit(opt.initStep,opt.Lip);

        if(opt.outLevel>=1) out.debug={}; end
        if(opt.adaptiveStep) cumu=0; end

        while(true)

            if(itr >= maxItr || (convThresh>2 && itr>=opt.minItr))
                break;
            end

            itr=itr+1;
            %if(mod(itr,100)==1 && itr>100) save('snapshotFST.mat'); end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  start of one PNPG step  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%

            numLineSearch=0; incStep=false;
            if(opt.adaptiveStep)
                if(cumu>=opt.cumuTol)
                    % adaptively increase the step size
                    t=t*opt.stepIncre;
                    cumu=0;
                    incStep=true;
                end
            %else
                %newTheta=(1+sqrt(1+4*B*theta^2))/2;
                if(itr==1)
                    newTheta=1;
                else
                    newTheta=1/gamma+sqrt(b+theta^2);
                end
                pbar=p+(theta-1)/newTheta*(p-preP);
                qbar=q+(theta-1)/newTheta*(q-preQ);
                [oldCost,gradp,gradq] = dualFunc(pbar,qbar);
            end

            % start of line Search
            while(true)
                numLineSearch = numLineSearch+1;

            %   if(opt.adaptiveStep)
            %       if(itr==1)
            %           newTheta=1;
            %       else
            %           B=t/preT;
            %           newTheta=1/gamma+sqrt(b+B*theta^2);
            %           %newTheta=(1+sqrt(1+4*B*theta^2))/2;
            %       end
            %       pbar=p+(theta-1)/newTheta*(p-preP);
            %       qbar=q+(theta-1)/newTheta*(q-preQ);
            %       [oldCost,gradp,gradq] = dualFunc(pbar,qbar);
            %   end

                [newP,newQ]=proxOp(pbar-gradp/t, qbar-gradq/t);

                newCost=dualFunc(newP,newQ);
                if((newCost-oldCost)<=...
                        innerProd(newP-pbar,gradp)+innerProd(newQ-qbar,gradq)...
                        +t*(sqrNorm(newP-pbar)+sqrNorm(newQ-qbar))/2)
                    if(itr<=opt.preSteps && opt.adaptiveStep && goodStep)
                        cumu=opt.cumuTol;
                    end
                    break;
                else
                    if(numLineSearch<=opt.maxLineSearch)
                        t=t/opt.stepShrnk; goodStep=false;
                        % Penalize if there is a step size increase just now
                        if(incStep)
                            incStep=false;
                            if(opt.incCumuTol)
                                opt.cumuTol=opt.cumuTol+4;
                            end
                        end
                    else  % don't know what to do, mark on debug and break
                        debug.appendLog('_FalseMM');
                        break;
                    end
                end
            end

            % using eps reduces numerical issue around the point of convergence
            if(newCost>cost && restart())
                itr=itr-1; continue;
            else
                difX = relativeDif([p(:);q(:)],[newP(:);newQ(:)]);
                preP = p; preQ = q;
                p = newP; q = newQ;
                theta = newTheta;
                preCost=cost;
                cost = newCost;
            end

            if(opt.adaptiveStep)
                preT=t;
                if(numLineSearch==1)
                    cumu=cumu+1;
                else
                    cumu=0;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %  end of one PNPG step  %
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            if(itr>1 && difX<=thresh )
                convThresh=convThresh+1;
            end

            if(opt.outLevel<1 && opt.debugLevel<2)
                continue;
            end

            difCost=abs(cost-preCost)/max(1,abs(cost));
            if(opt.outLevel>=1)
                out.time(itr)=toc(tStart);
                out.cost(itr)=cost;
                out.difX(itr)=difX;
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
                debug.print(2,sprintf(' %14.12g',cost));
                debug.print(2,sprintf(' %12g %4d',difX,numLineSearch));
                debug.print(2,sprintf(' %12g', difCost));
                debug.clear_print(2);
                if(mod(itr,opt.verbose)==0) debug.println(2); end
            end
        end
        Qp=a-u*(Psi_v(p)+Psi_h(q)); % is real
        x=prj_C(Qp);   % is real
        pOut={p,q};
        out.opt = opt; out.date=datestr(now);
        if(opt.outLevel>=2)
            out.gradp=gradp;
            out.gradq=gradq;
        end
        if(debug.level(1))
            fprintf('\nCPU Time: %g, objective=%g\n',toc(tStart),cost);
        end

        function res=restart()
            % if has monmentum term, restart
            res=((pNorm(pbar-p,1)+pNorm(qbar-q,1))~=0);
            if(res)
                theta=1;
                debug.appendLog('_Restart');
                debug.printWithoutDel(2,'\t restart');
            end
        end

        function t=stepSizeInit(select,Lip,delta)
            switch (lower(select))
                case {'bb','hessian'}   % use BB method to guess the initial stepSize
                    if(~exist('delta','var')) delta=1e-5; end
                    [~,gradp1,gradq1] = dualFunc(p,q);
                    tt=sqrt(sqrNorm(gradp1)+sqrNorm(gradq1));
                    tempp= delta*gradp1/tt;
                    tempq= delta*gradq1/tt;
                    [~,gradp2,gradq2] = dualFunc(p-tempp,q-tempq);
                    t = abs(innerProd(gradp1-gradp2,tempp)+innerProd(gradq1-gradq2,tempq))/delta^2;
                case 'fixed'
                    t = Lip;
                otherwise
                    error('unkown selection for initial step');
            end
            if(isnan(t) || t<=0)
                error('\n PNPG is having a negative or NaN step size, do nothing and return!!\n');
            end
        end

    end
    function [D,iter,pOut,out]=denoise(Xobs,lambda,epsilon,MAXITER,pInit)
        % modified by renliang gu
        %
        %This function implements the FISTA method for TV-based denoising problems
        %
        % Based on the paper
        % Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
        % Total Variation Image Denoising and Deblurring Problems"
        % -----------------------------------------------------------------------
        % Copyright (2008): Amir Beck and Marc Teboulle
        % 
        % FISTA is distributed under the terms of 
        % the GNU General Public License 2.0.
        % 
        % Permission to use, copy, modify, and distribute this software for
        % any purpose without fee is hereby granted, provided that this entire
        % notice is included in all copies of any software which is or includes
        % a copy or modification of this software and in all copies of the
        % supporting documentation for such software.
        % This software is being provided "as is", without any express or
        % implied warranty.  In particular, the authors do not make any
        % representation or warranty of any kind concerning the merchantability
        % of this software or its fitness for any particular purpose."
        % ----------------------------------------------------------------------
        % INPUT
        % Xobs ..............................an observed noisy image.
        % lambda ........................ parameter
        %  
        % OUTPUT
        % D ........................... The solution of the problem 
        %                                            min{||X-Xobs||^2+2*lambda*TV(X
        %                                            ) : l <= X_{i,j} <=u} 
        % iter .............................  Number of iterations required to get
        %                                            an optimal solution (up to a tolerance)
        % fun_all ......................   An array containing all the function
        %                                             values obtained during the
        %                                             iterations

        %Define the Projection onto the box

        if (~exist('MAXITER','var'))
            MAXITER=100;
        end
        if (~exist('epsilon','var'))
            epsilon=1e-4;
        end

        [m,n]=size(Xobs);
        D=zeros(m,n);
        if(exist('pInit','var') && ~isempty(pInit))
            P1=pInit{1}; P2=pInit{2};
        else
            [P1,P2]=Ltrans(D);
        end

        switch tv
            case 'iso'
                A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
                A=sqrt(max(A,1));
                P1=P1./A(1:m-1,:);
                P2=P2./A(:,1:n-1);
            case {'l1', '1d'}
                P1=max(min(P1,1),-1);
                P2=max(min(P2,1),-1);
            otherwise
                error('unknown type of total variation. should be iso or l1');
        end
        R1=P1; R2=P2;

        % debug=true;
        f=@(x,u) 0.5*norm(x-Xobs,'fro')^2+u*tlv(x,tv);

        tk=1;
        tkp1=1;
        count=0;
        i=0;

        fval=inf;
        fun_all=[];
        if(prnt)
            fprintf('***********************************\n');
            fprintf('*Solving with FGP/FISTA**\n');
            fprintf('***********************************\n');
            fprintf('#iteration  function-value  relative-difference\n');
            fprintf('---------------------------------------------------------------------------------------\n');
        end
        while((i<MAXITER)&&(count<5))
            fold=fval;
            %%%%%%%%%
            % updating the iteration counter
            i=i+1;
            %%%%%%%%%
            % Storing the old value of the current solution
            Dold=D;
            %%%%%%%%%%
            %Computing the gradient of the objective function
            Pold1=P1;
            Pold2=P2;
            tk=tkp1;
            C=Xobs-lambda*Lforward(R1,R2);
            D=prj_C(C);
            [Q1,Q2]=Ltrans(D);

            %%%%%%%%%%
            % Taking a step towards minus of the gradient
            P1=R1+Q1/(8*lambda);
            P2=R2+Q2/(8*lambda);

            %%%%%%%%%%
            % Peforming the projection step
            switch tv
                case 'iso'
                    A=[P1;zeros(1,n)].^2+[P2,zeros(m,1)].^2;
                    A=sqrt(max(A,1));
                    P1=P1./A(1:m-1,:);
                    P2=P2./A(:,1:n-1);
                case {'l1','1d'}
                    P1=max(min(P1,1),-1);
                    P2=max(min(P2,1),-1);
                    %     P1=P1./(max(abs(P1),1));
                    %     P2=P2./(max(abs(P2),1));
                otherwise
                    error('unknown type of total variation. should be iso or l1');
            end

            %%%%%%%%%%
            %Updating R and t
            tkp1=(1+sqrt(1+4*tk^2))/2;

            R1=P1+(tk-1)*(P1-Pold1)/tkp1;
            R2=P2+(tk-1)*(P2-Pold2)/tkp1;

            re=norm(D-Dold,'fro')/norm(D,'fro');
            if (re<epsilon)
                count=count+1;
            else
                count=0;
            end

            fval=(-norm(C-D,'fro')^2+norm(C,'fro')^2)/2;
            if (fval>fold)
                tkp1=1;
            end

            if(nargout==4)
                fun_all=[fun_all;fval];
                if(prnt)
                    fprintf('%7d %10.10g %10.10g  %19g',i,fval,norm(D-Dold,'fro')/norm(D,'fro'),f(D,lambda));
                    if (fval>fold) fprintf('  *\n'); else fprintf('   \n'); end
                end
            end
        end
        iter=i;
        pOut={P1, P2};
        out.cost=fun_all;
    end

end

% If edit the following, update TV.[A,B]t\=
function x = Psi_v(p)
    [I,J]=size(p);
    I=I+1;
    x=[p; zeros(1,J)];
    x(2:I,:)=x(2:I,:)-p(1:I-1,:);
end
function p = Psi_vt(x)
    [I,~]=size(x);
    p=x(1:I-1,:)-x(2:I,:);
end
function x = Psi_h(q)
    [I,J]=size(q);
    J=J+1;
    x=[q zeros(I,1)];
    x(:,2:J)=x(:,2:J)-x(:,1:J-1);
end
function q = Psi_ht(x)
    [~,J]=size(x);
    q=x(:,1:J-1)-x(:,2:J);
end
function [p,q]=isoPrj(p,q)
    [I,J]=size(p);
    I=I+1;
    %mag=sqrt(max(1,[p.^2;zeros(1,J)]+[q.^2, zeros(I,1)]));
    mag=sqrt(max(1,p(:,1:J-1).^2+q(1:I-1,:).^2));
    p(:,1:J-1)=p(:,1:J-1)./mag; p(:,J)=min(max(p(:,J),-1),1);
    q(1:I-1,:)=q(1:I-1,:)./mag; q(I,:)=min(max(q(I,:),-1),1);
end
function [p,q]=l1Prj(p,q)
    p=min(max(p,-1),1);
    q=min(max(q,-1),1);
end

