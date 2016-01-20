function [alpha,pppp] = admm(Psi,Psit,a,u,relativeTol,maxItr,isInDebugMode,init)
    %
    % solve 0.5*||α-a||_2 + I(α≥0) + u*||Psit(α)||_1
    %
    % author: Renliang Gu (gurenliang@gmail.com)
    %
    if((~exist('relativeTol','var')) || isempty(relativeTol)) relativeTol=1e-6; end
    if((~exist('maxItr','var')) || isempty(maxItr)) maxItr=1e3;  end
    if((~exist('init','var')) || isempty(init)) init=a;  end
    if((~exist('isInDebugMode','var')) || isempty(isInDebugMode)) isInDebugMode=false;  end
    % this makes sure the convergence criteria is nontrival
    relativeTol=min(1e-3,relativeTol);
    nu=0; rho=1; cnt=0; preS=Psit(init); s=preS;

    pppp=0;
    while(true)
        pppp=pppp+1;
        cnt= cnt + 1;

        alpha = max((a+rho*Psi(s+nu))/(1+rho),0);
        Psit_alpha=Psit(alpha);
        s = Utils.softThresh(Psit_alpha-nu,u/rho);
        nu=nu+s-Psit_alpha;

        difS=pNorm(s-preS); preS=s;
        residual = pNorm(s-Psit_alpha);
        sNorm = max(pNorm(s),eps);

        if(isInDebugMode)
            cost(pppp)=0.5*sqrNorm(max(Psi(s),0)-a)+u*pNorm(Psit(max(Psi(s),0)),1);
            cost1(pppp)=0.5*sqrNorm(alpha-a)+u*pNorm(Psit_alpha,1);
            cost2(pppp)=0.5*sqrNorm(alpha-a)+u*pNorm(s,1);

            rhoRec(pppp)=rho;

            if(pppp>1)
                difAlpha = pNorm(preAlpha-alpha);
                difAlphaRec(pppp-1)=difAlpha/sNorm;
                difSRec(pppp-1)=difS/sNorm;
                residualRec(pppp-1)=residual/sNorm;
            end

            preAlpha=alpha;
        end

        if(pppp>maxItr) break; end
        if(difS<=relativeTol*sNorm && residual<=relativeTol*sNorm) break; end
        if(cnt>10) % prevent excessive back and forth adjusting
            if(difS>10*residual)
                rho = rho/2 ; nu=nu*2; cnt=0;
            elseif(difS<residual/10)
                rho = rho*2 ; nu=nu/2; cnt=0;
            end
        end
    end 
    alpha = max((a+rho*Psi(s+nu))/(1+rho),0);

    if(isInDebugMode)
        costRef=0.5*sqrNorm(max(init,0)-a)+u*pNorm(Psit(max(init,0)),1);
        figure;
        semilogy(rhoRec,'r'); title('rho');
        figure;
        semilogy(difAlphaRec,'r'); hold on;
        semilogy(difSRec,'g');
        semilogy(residualRec,'b');
        legend('difAlpha','difS','residual');
        title('covnergence criteria');
        figure;
        semilogy(cost1-min([cost,cost1,cost2]),'r-.'); hold on;
        semilogy(cost -min([cost,cost1,cost2]),'g-'); hold on;
        semilogy(cost2 -min([cost,cost1,cost2]),'b-'); hold on;
        title('admm centered obj');
        legend('alpha','s','alpha,s');
        figure;
        semilogy(cost1,'r'); hold on;
        semilogy(cost,'g'); hold on;
        semilogy(ones(size(cost))*costRef,'k');
        title('admm obj');
        legend('alpha','s','ref');
    end
    % end of the ADMM inside the NPG
end
