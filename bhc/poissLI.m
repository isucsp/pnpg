function [f,g,h] = poissLI(Imea,A,Ie)
    Ir=A*Ie; Err=log(Ir./Imea);
    f=sum(Ir-Imea)-Imea'*Err;
    if(nargout>1) g=A'*(1-Imea./Ir); end
    if(nargout>2)
        weight=Imea./(Ir.^2);
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

