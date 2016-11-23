classdef ProximalADMM < Proximal
    %
    % solve 0.5*||α-a||_2^2 + I(α≥0) + u*||Psit(α)||_1
    %
    % author: Renliang Gu (gurenliang@gmail.com)
    %
    properties
        opt
        prj_C=@(x)x;
        Psi
        Psit
        out % put this property to be private later
        isInDebugMode
    end
    properties (Access=private)
        initNu
        zeroInitSize
        temp
        nu
        s
    end
    properties (Dependent)
        tvType
    end

    methods
        function obj = ProximalADMM(Psi,Psit,prj_C)
            ??? test psi and Psit and incases that they are matrix
            obj.Psi=Psi;
            obj.Psit=Psit;
            if(exist('prj_C','var'))
                obj.prj_C=prj_C;
            else
                obj.prj_C=@(x) x;
            end
            obj.iterative=true;
            obj.preConfig();

            obj.initNu=0;

            ???
            obj.zeroInitSize=@(x) size(x);
            obj.regFunc = @(X) tlv(X,tvType);
        end
        % don't forget to set default value for maxItr, thresh
        function preConfig(obj)
            obj.thresh=1e-6; obj.maxItr=1e3; obj.isInDebugMode=false;
        end
        function setInit(obj)
            obj.init=obj.s;  obj.initNu=obj.nu;
        end
        function set.thresh(obj,thresh)
            % this makes sure the convergence criteria is nontrival
            obj.thresh=min(1e-3,obj.thresh);
        end
        function alpha=denoise(obj,a,u)
            obj.config(a,u);
            strlen=0;

            obj.x=obj.prj_C(a-u*obj.temp);
            obj.steps=obj.out.p;
            x=obj.x;

            % scale the input to prevent numerical problem
            scale=pNorm(a);
            if(scale==0)
                alpha=zeros(size(a)); obj.steps=0;
                return;
            end
            a=a/scale;  u=u/scale;
            nu=obj.initNu/scale; preS=obj.init/scale; s=preS;

            rho=1; cnt=0;

            pppp=0;
            while(true)
                pppp=pppp+1;
                cnt= cnt + 1;

                alpha = obj.prj_C((a+rho*obj.Psi(s+nu))/(1+rho));
                Psit_alpha=obj.Psit(alpha);
                s = Utils.softThresh(Psit_alpha-nu,u/rho);
                nu=nu+s-Psit_alpha;

                difS=pNorm(s-preS); preS=s;
                residual = pNorm(s-Psit_alpha);

                if(obj.isInDebugMode)
                    cost=0.5*sqrNorm(alpha-a)+u*pNorm(Psit_alpha,1);
                    gap=rho*nu'*(s-Psit_alpha);

                    str=sprintf('itr=%d, cost=%g pRes=%g dRes=%g gap=%g rho=%g       ',pppp,...
                        cost,residual,difS,gap,rho);
                    if(strlen==0 || (mod(pppp-1,100)==0 || (pppp<=100 && mod(pppp-1,10)==0) || pppp-1<10))
                        fprintf('\n%s',str);
                    else
                        fprintf([repmat('\b',1,strlen) '%s'],str);
                    end
                    strlen = length(str);
                end

                if(pppp>obj.maxItr) break; end
                if(difS<=obj.thresh && residual<=obj.thresh) break; end
                if(cnt>20) % prevent excessive back and forth adjusting
                    if(difS>10*residual)
                        rho = rho/2 ; nu=nu*2; cnt=0;
                    elseif(difS<residual/10)
                        rho = rho*2 ; nu=nu/2; cnt=0;
                    end
                end
            end 
            alpha = obj.prj_C((a+rho*obj.Psi(s+nu))/(1+rho))*scale;
            obj.x=alpha;
            obj.nu=nu*rho*scale;
            obj.s =s*scale;
            obj.Psit_alpha=Psit_alpha;
        end
        function y=getPenalty(obj,X)
            if(~exist('X','var'))
                y=pNorm(obj.Psit_alpha,1);
            else
                y=pNorm(obj.Psit(X),1);
            end
        end
        %     ??? check complex variable of solver
        function config(obj,a,u)
            obj.opt.thresh=obj.thresh; obj.opt.maxItr=obj.maxItr;
            if(any(size(obj.init)-obj.zeroInitSize(a)~=0))
                obj.init=zeros(obj.zeroInitSize(a));
            end
            obj.opt.NLL=@(p) dualTvFunc(p,a,obj.Psi,obj.Psit,u,obj.prj_C);
            function [f,g,h] = dualTvFunc(p,a,Psi,Psit,u,prj_C)
                % minimize 0.5*||a-u*real(Psi(p))||_2^2-0.5*||X-a+u*real(Psi(p))||_2^2
                % subject to ||p||_infty <= 1
                % where X=prj_C( a-u*real(Psi(p)) ), p and Psi may be complex.
                Qp=a-u*Psi(p); % is real
                x=prj_C(Qp);   % is real
                f=(sqrNorm(Qp)-sqrNorm(x-Qp))/2;
                if(nargout>1)
                    g=-u*Psit(x);
                    if(nargout>2)
                        h=[];
                    end
                end
            end
        end
    end
end

end

% if(obj.isInDebugMode)
%     rhoRec(pppp)=rho;
% 
%     if(pppp>1)
%         difAlpha = pNorm(preAlpha-alpha);
%         difAlphaRec(pppp-1)=difAlpha/sNorm;
%         difSRec(pppp-1)=difS/sNorm;
%         residualRec(pppp-1)=residual/sNorm;
%     end
% 
%     preAlpha=alpha;
% 
%     costRef=0.5*sqrNorm(max(init,0)-a)+u*pNorm(Psit(max(init,0)),1);
%     figure;
%     semilogy(rhoRec,'r'); title('rho');
%     figure;
%     semilogy(difAlphaRec,'r'); hold on;
%     semilogy(difSRec,'g');
%     semilogy(residualRec,'b');
%     legend('difAlpha','difS','residual');
%     title('covnergence criteria');
%     figure;
%     semilogy(cost1-min([cost,cost1,cost2]),'r-.'); hold on;
%     semilogy(cost -min([cost,cost1,cost2]),'g-'); hold on;
%     semilogy(cost2 -min([cost,cost1,cost2]),'b-'); hold on;
%     title('admm centered obj');
%     legend('alpha','s','alpha,s');
%     figure;
%     semilogy(cost1,'r'); hold on;
%     semilogy(cost,'g'); hold on;
%     semilogy(ones(size(cost))*costRef,'k');
%     title('admm obj');
%     legend('alpha','s','ref');
% end
