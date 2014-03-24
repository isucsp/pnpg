function [f,g,h] = gaussLAlpha(Imea,Ie,alpha,Phi,Phit,polyIout,obj)
    for i=1:size(alpha,2) R(:,i)=Phi(alpha(:,i)); end
    if(nargout>=3)
        [BLI,sBLI,ssBLI] = polyIout(R,Ie); %A=exp(-R*mu');
    elseif(nargout>=2)
        [BLI,sBLI] = polyIout(R,Ie); %A=exp(-R*mu');
    else
        BLI = polyIout(R,Ie); %A=exp(-R*mu');
    end
    Err=log(BLI./Imea);
    f=Err'*Err/2;
    if(nargin>=7)
        obj.zmf=[min(Err(:)); max(Err(:))]; % lb and ub of z-f(theta)
    end
    if(nargout>=2)
        g=-Phit(Err.*sBLI./BLI);
        if(nargout>=3)
            weight=(1-Err).*((sBLI./BLI).^2)+Err.*ssBLI./BLI;
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

