function [f,g,h] = poissLAlpha(Imea,Ie,alpha,Phi,Phit,polyIout)
    for i=1:size(alpha,2) R(:,i)=Phi(alpha(:,i)); end
    if(nargout>=3)
        [BLI,sBLI,ssBLI] = polyIout(R,Ie); %A=exp(-R*mu');
    elseif(nargout>=2)
        [BLI,sBLI] = polyIout(R,Ie); %A=exp(-R*mu');
    else
        BLI = polyIout(R,Ie); %A=exp(-R*mu');
    end
    Err=log(BLI./Imea);
    f=sum(BLI-Imea)-Imea'*Err;
    if(nargout>=2)
        g=Phit((1-Imea./BLI).*sBLI);
        if(nargout>=3)
            weight=((sBLI./BLI).^2).*Imea+ssBLI.*(1-Imea./BLI);
            h=@(x,opt) hessian(Phi, Phit, weight, x,opt);
        end
    end
end
function h = hessian(Phi, Phit, weight, x, opt)
    y = Phi(x);
    if(opt==1)
        h = Phit(weight.*y);
    else
        h= weight'*(y.*y);
    end
end

