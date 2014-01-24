
function [f,zmf,g,h] = fIe(Imea,A,Ie)
    % Err= z-f(theta)
    Ir=A*Ie; Err=log(Ir./Imea); f=Err'*Err;
    zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    if(nargout>2) g=2*A'*(Err./Ir); end
    if(nargout>3)
        %Err=Err*0;
        h=2*A'*(repmat((1-Err)./(Ir.^2),1,length(Ie)).*A);
    end
end

