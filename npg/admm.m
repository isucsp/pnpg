function [x,itr,p] = admm(Psi,Psit,a,u,relativeTol,maxItr,isInDebugMode,pInit)
    %
    % solve 0.5*||α-a||_2^2 + I(α≥0) + u*||Psit(α)||_1
    %
    % author: Renliang Gu (gurenliang@gmail.com)
    %
    if((~exist('relativeTol','var')) || isempty(relativeTol)) relativeTol=1e-6; end
    if((~exist('maxItr','var')) || isempty(maxItr)) maxItr=1e3;  end
    if((~exist('isInDebugMode','var')) || isempty(isInDebugMode)) isInDebugMode=false;  end
    if((~exist('pInit','var')) || isempty(pInit) || ~iscell(pInit))
        temp=size(Psit(a));
        pInit={zeros(temp), zeros(temp), 1};
    end
    % this makes sure the convergence criteria is nontrival
    relativeTol=min(1e-3,relativeTol); scale=1; strlen=0;

    % scale the input to prevent numerical problem
    scale=pNorm(a,2); if(scale==0) x=zeros(size(a)); itr=0; p=pInit; return; end
    s=pInit{1}; nu=pInit{2}; rho=pInit{3};

    a=a/scale;  u=u/scale; s=s/scale; nu=nu/scale;

    itr=0; cnt=0; preS=s;
    while(true)
        itr=itr+1;
        cnt= cnt + 1;

        x = max((a+rho*Psi(s+nu))/(1+rho),0);
        Psit_x=Psit(x);
        s = Utils.softThresh(Psit_x-nu,u/rho);
        nu=nu+s-Psit_x;

        difS=pNorm(s-preS,2); preS=s;
        residual = pNorm(s-Psit_x,2);

        if(isInDebugMode)
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
    x = max(scale*(a+rho*Psi(s+nu))/(1+rho),0);
    p = {s,nu,rho};
    % end of the ADMM inside the NPG
end



% if(isInDebugMode)
%     rhoRec(itr)=rho;
% 
%     if(itr>1)
%         difAlpha = pNorm(preAlpha-x,2);
%         difAlphaRec(itr-1)=difAlpha/sNorm;
%         difSRec(itr-1)=difS/sNorm;
%         residualRec(itr-1)=residual/sNorm;
%     end
% 
%     preAlpha=x;
% 
%     costRef=0.5*sqrNorm(max(pInit,0)-a)+u*pNorm(Psit(max(pInit,0)),1);
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
%     legend('x','s','x,s');
%     figure;
%     semilogy(cost1,'r'); hold on;
%     semilogy(cost,'g'); hold on;
%     semilogy(ones(size(cost))*costRef,'k');
%     title('admm obj');
%     legend('x','s','ref');
% end
