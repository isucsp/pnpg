function [BL,sBL,ssBL] = disIout(kappa, s, I)

    if(nargin==0)       % start to test
        kappa = logspace(log10(0.03),log10(30),100);
        s = -0.1:0.001:10;
        [BL, sBL, ssBL] = disIout(kappa,s);
        return;
    end

    kappa = kappa(:); s = s(:);
    if(nargin>=3)
        BL = zeros(length(s),1);
    else
        BL = zeros(length(s),length(kappa)-1);
    end
    if(nargout>1) sBL = BL; end
    if(nargout>2) ssBL = BL; end

    for i=1:length(kappa)
        temp = exp(-kappa(i)*s);
        if(nargin>=3) BL = BL + temp*I(i);
        else BL(:,i) = temp; end

        if(nargout>1)
            temp = temp*kappa(i);
            if(nargin>=3) sBL = sBL + temp*I(i);
            else sBL(:,i) = temp; end
        end
        if(nargout>2)
            temp = temp*kappa(i);
            if(nargin>=3) ssBL = ssBL + temp*I(i);
            else ssBL(:,i) = temp; end
        end
    end
end

