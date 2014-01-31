
function [f,g,h] = gaussLI(Imea,A,Ie)
    % Err= z-f(theta)
    Ir=A*Ie; Err=log(Ir./Imea); f=Err'*Err;
    if(nargout>1) g=2*A'*(Err./Ir); end
    if(nargout>2)
        %Err=Err*0;
        h=2*A'*(repmat((1-Err)./(Ir.^2),1,length(Ie)).*A);
    end
    %if(nargout>3)
    %    zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    %end
end

