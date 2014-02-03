function [f,zmf,g,weight] = gaussLAlpha(Imea,Ie,alpha,kappa,Phi,Phit,polyIout)
    for i=1:size(alpha,2) R(:,i)=Phi(alpha(:,i)); end
    if(nargout>3)
        [BI,temp,ssBLI] = polyIout(R,Ie); %A=exp(-R*mu');
    elseif(nargout>2)
        [BI,temp] = polyIout(R,Ie); %A=exp(-R*mu');
    else
        BI = polyIout(R,Ie); %A=exp(-R*mu');
    end
    Err=log(BI./Imea); f=Err'*Err;
    zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    if(nargout>=3)
        g=-2*Phit(Err.*temp./BI);
        if(nargout>=4)
            weight=(1-Err).*((temp./BI).^2)+Err.*ssBLI./BI;
        end
    end
end

