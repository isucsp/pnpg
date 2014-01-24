
function [f,zmf,g,weight] = fAlpha(Imea,Ie,R,mu,Phi,Phit,tdphi)
    if(nargout>3)
        [A,temp,ssBLI] = polyOutput('b0',kappa,R,Ie); %A=exp(-R*mu');
    else if(nargout>2)
        [A,temp] = polyOutput('b0',kappa,R,Ie); %A=exp(-R*mu');
    else
        A = polyOutput('b0',kappa,R,Ie); %A=exp(-R*mu');
    end
    Ir=A*Ie; Err=log(Ir./Imea); f=Err'*Err;
    zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    if(nargout>2)
        g=-2*Phit(Err.*temp./Ir);
        weight=(1-Err).*((temp./Ir).^2)+Err.*ssBLI./Ir;
        %if(min(weight)<=0)
        %    fprintf([repmat(' ',1,80) repmat('\n',1,1)]);
        %    fprintf(['fAlpha[warning]: obj is non-convex over alpha,'...
        %        'minimum=%g\n'],min(weight));
        %    fprintf([repmat(' ',4,80) repmat('\n',4,1)]);
        %end
        %temp=Phi(g+tdphi);
        %h=2*temp'*(weight.*temp);
    end
end

