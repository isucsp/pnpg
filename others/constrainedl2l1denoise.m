function [x,iter] = constrainedl2l1denoise(y,W,WT,tau,mu,miniter,maxiter,...
    stopcriterion,tolerance,isInDebugMode)
% For now, the tolerance is the relative duality gap, and is the only
% convergence criterion implemented
% Also, in the future it would be good to output the number of nonzeros
% in theta, even though we output a solution in x.
% May be worthwhile to clean this up and keep as a separate function,
% however this would entail coding such that it checks the inputs and
% outputs and such...

verbose = 0; % For now, do not output any status
if((~exist('isInDebugMode','var')) || isempty(isInDebugMode)) isInDebugMode=false;  end

gamma = zeros(size(y));
lambda = gamma;
WTlambda = WT(lambda);
y = WT(y);
iter = 1;
converged = 0;

while (iter <= miniter) || ((iter <= maxiter) && not(converged))
    %disp(['Subiter = ',num2str(iter)])
    gamma       = min( max( -tau, -y - WTlambda), tau);
    lambda      = max( -W(y + gamma) - mu, 0);
    WTlambda    = WT(lambda);
    theta       = y + gamma + WTlambda;

    if(isInDebugMode)
        temp=abs(W(theta));
        cost(iter)=0.5*sqrNorm(temp-W(y))+tau*pNorm(WT(temp),1);
        if(iter>1)
            difTheta(iter) = pNorm(preTheta-theta)/pNorm(theta);
        end
    end

    % Check for convergence
    if iter >= miniter % no need to check if miniter not reached
        switch stopcriterion
            case 0
                % Just exhaust maxiter
                converged = 0;
            case 1
                primal_obj = sum( (theta(:)-y(:)).^2)./2 + tau.*sum(abs(theta(:)));
                dual_obj = -sum( theta(:).^2)./2 + sum(y(:).^2)./2-mu.*sum(lambda(:));
                % Need this for what's in verbose:
                % duality_gap = primal_obj - dual_obj;
                rel_duality_gap = abs(primal_obj-dual_obj)/max(-primal_obj,dual_obj);
                if verbose
                    % display some stuff
                    % fprintf('l1Den: It=%4d, PObj=%13.5e, DObj=%13.5e, DGap=%13.5e, RDGap=%13.5e\n', iter,primal_obj,dual_obj,duality_gap,rel_duality_gap)
                end
                if (rel_duality_gap <= tolerance) || isinf(rel_duality_gap)
                    converged = 1;
                end
            case 2
                if(iter>1 && pNorm(preTheta-theta)<=tolerance*max(1,pNorm(theta)))
                    converged = 1;
                end
                preTheta=theta;
        end
    end
    iter = iter + 1;
end    
x = abs(W(theta)); %note, sometimes W returns small negative values

if(isInDebugMode)
    figure;
    semilogy(cost-min(cost),'g.'); hold on;
    title('spiral centered cost');
    figure;
    semilogy(cost,'b'); hold on;
    title('spiral cost');
    figure;
    semilogy(difTheta,'b'); hold on;
    title('spiral diftheta');
end
end

