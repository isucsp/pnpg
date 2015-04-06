
function [f,g,h] = gaussLI(Imea,A,Ie)
    % Err= z-f(theta)
    Ir=A*Ie; Err=log(Ir./Imea); f=Err'*Err/2;
    if(nargout>1) g=A'*(Err./Ir); end
    if(nargout>2)
        %Err=Err*0;
        %h=A'*(repmat((1-Err)./(Ir.^2),1,length(Ie)).*A);
        weight=(1-Err)./(Ir.^2);
        h = @(x,opt) hessian(x,opt);
    end
    %if(nargout>3)
    %    zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    %end
    function h = hessian(x,opt)
        y = A*x;
        if(opt==1)
            h = A'*(weight.*y);
        elseif(opt==2) h = weight'*(y.*y);
        else
            n=length(weight);
            h=y'*sparse(1:n,1:n,weight,n,n)*y;
        end
    end
end

