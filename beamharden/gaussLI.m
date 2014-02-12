function [f,g,h] = gaussLI(Imea,A,Ie)
    % Err= z-f(theta)
    Ir=A*Ie; Err=log(Ir./Imea); f=Err'*Err/2;
    if(nargout>1) g=A'*(Err./Ir); end
    if(nargout>2)
        %Err=Err*0;
        %h=2*A'*(repmat((1-Err)./(Ir.^2),1,length(Ie)).*A);
        h = @(x,opt) hessian(A,(1-Err)./(Ir.^2),x,opt);
    end
    %if(nargout>3)
    %    zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    %end
end

function h = hessian(A,weight,x,opt)
    y = A*x; h=(repmat(weight,1,size(x,2)).*y);
    if(opt==1)
        h = A'*h;
    else
        h = y'*h;
    end
end

