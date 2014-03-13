classdef Utils < handle
    methods(Static)
        function y = softThresh(x,thresh)
            idx = abs(x)<thresh;
            y = zeros(size(x));
            y(x>0) = x(x>0)-thresh;
            y(x<0) = x(x<0)+thresh;
            y(idx) = 0;
        end
        function [f,g,h] = augLag(x,z)
            g=x-z;
            f=norm(x-z)^2/2;
            if(nargout>=3)
                h = @(xx,opt) hessian(xx,opt);
            end
            function hh = hessian(x,opt)
                if(opt==1) hh = x;
                else hh = x'*x;
                end
            end
        end
        function [f,g,h] = nonnegLogBarrier(alpha)
            %if(any(alpha(:)<=0)) f=eps^-1; alpha(alpha<=0)=eps;
            %else f=-sum(log(alpha(:)));
            %end
            %if(nargout>1) g = -1./alpha; h=1./alpha.^2; end
            f=-sum(log(alpha(:)));
            if(nargout>1)
                g=-1./alpha;
                h=1./(alpha.^2);
            end
        end
        function [f,g,h] = nonnegPen(alpha)
            temp=(alpha<0);
            f=alpha(temp)'*alpha(temp);
            if(nargout>=2)
                g=zeros(size(alpha));
                g(temp)=2*alpha(temp);
                if(nargout>=3)
                    h = @(x,opt) hessian(x,opt);
                end
            end
            function hh = hessian(x,opt)
                if(opt==1)
                    hh = zeros(size(x));
                    hh(temp,:) = x(temp,:)*2;
                else
                    y = x(temp,:);
                    hh = y'*y*2;
                end
            end
        end

        function [f,g,h]=barrierIe(Ie)
            %if(any(Ie)<=0)
            %    Ie(Ie<=0)=eps; f=eps^-1;
            %    if(1-sum(Ie)<=0) Ie=Ie*(1-eps)/sum(Ie); end
            %else
            %    if(1-sum(Ie)<=0) Ie=Ie*(1-eps)/sum(Ie); f=eps^-1;
            %    else f=-sum(log(Ie))-log(1-sum(Ie)); end
            %end
            f=-sum(log(Ie))-log(1-sum(Ie));
            if(nargout>1)
                g=1/(1-sum(Ie))-1./Ie;
                h=1/(1-sum(Ie))^2+diag(1./(Ie.^2));
            end
        end

        function [f,g,h] = lustigL1(alpha,xi,Psi,Psit)
            s=Psit(alpha);
            sqrtSSqrMu=sqrt(s.^2+xi);
            f=sum(sqrtSSqrMu);
            if(nargout>=2)
                g=Psi(s./sqrtSSqrMu);
                if(nargout>=3)
                    h = @(x,opt) hessian(xi./(sqrtSSqrMu.^3),x,opt);
                end
            end
            function hh = hessian(weight,x,opt)
                y = Psit(x);
                if(opt==1)
                    hh = Psi(weight.*y);
                else
                    hh = y'*(weight.*y);
                end
            end
        end

        function [f,g,h] = huber(alpha,mu,Psi,Psit)
            s=Psit(alpha);
            idx = abs(s)<mu;
            temp = abs(s)-mu/2;
            temp(idx) = s(idx).^2/2/mu;
            f = sum(temp);
            if(nargout>=2)
                temp = ones(size(s));
                temp(idx) = s(idx)/mu;
                g=Psi(temp);
                if(nargout>=3)
                    h = @(x,opt) hessian(x,opt);
                end
            end
            function hh = hessian(x,opt)
                y = Psit(x);
                if(opt==1)
                    y(idx) = y(idx)/mu; y(~idx) = 0; hh=Psi(y);
                else
                    y = y(idx); hh = y'*y/mu;
                end
            end
        end

        function [f,g,h] = linearModel(alpha,Phi,Phit,y)
            PhiAlpha=Phi(alpha);
            f=norm(y-PhiAlpha)^2/2;
            if(nargout>=2)
                g=Phit(PhiAlpha-y);
                if(nargout>=3)
                    h=@(x,opt) hessian(x,opt);
                end
            end
            function hh = hessian(x,opt)
                yy = Phi(x);
                if(opt==1)
                    hh = Phit(yy);
                else
                    hh = yy'*yy;
                end
            end
        end

    end
end
