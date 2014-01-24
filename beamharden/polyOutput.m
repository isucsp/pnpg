function [BL,sBL,ssBL] = polyOutput(option, kappa, s, I)
    % kappa is a column vector, representing the nodes
    
    if(nargin==0)       % start to test
        kappa = logspace(log10(0.03),log10(30),100);
        s = -0.1:0.001:10;
        option = 'b0';
        nargout = 3;
    end
    
    kappa = kappa(:); s = s(:);
    idx = find(abs(s)<=eps);
    switch lower(option)
        case 'b0'       % B-0 spline with nodes be kappa
            if(nargin>=4)
                BL
            BL = zeros(length(s),length(kappa)-1);
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end
            expksR = exp(-kappa(1)*s);
            for i=1:length(kappa)-1
                expksL = expksR;
                expksR = exp(-kappa(i+1)*s);
                BL(:,i) = (expksL-expksR)./s;
                BL(idx,i) = kappa(i+1)-kappa(i);
                if(nargout>1)
                    sBL(:,i) = ((s*kappa(i)+1).*expksL-(s*kappa(i+1)+1).*expksR)./s./s;
                    sBL(idx,i) = (kappa(i+1)-kappa(i))*(kappa(i+1)+kappa(i))/2;
                end
                if(nargout>2)
                    ssBL(:,i) = ( ((kappa(i)*s+2)*kappa(i).*s+2).*expksL ...
                        -((kappa(i+1)*s+2)*kappa(i+1).*s + 2).*expksR )./s./s./s;
                    ssBL(idx,i) = (kappa(i+1)^3-kappa(i)^3)/3;
                end
            end
%        case 'b1'

    end
end

