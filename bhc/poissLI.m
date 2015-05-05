function [f,g,h] = poissLI(Imea,A,Ie)
    Ir=A*Ie(:);
    if(nargout>1)
        ImeaOverIr=Imea./Ir;
        Err=-log(ImeaOverIr);
        g=A'*(1-ImeaOverIr);
    else
        Err=log(Ir./Imea);
    end
    f=sum(Ir-Imea)-Imea'*Err;
    if(nargout>2)
        weight=ImeaOverIr./Ir;
        h = @(x,opt) hessian(x,opt);
    end
    function h = hessian(x,opt)
        y = A*x;
        if(opt==1)
            h = A'*(weight.*y);
        else
            h=[];
            for i=1:size(y,2)
                h(:,i)=y'*(y(:,i).*weight);
            end
        end
    end
end

