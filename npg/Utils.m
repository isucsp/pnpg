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
            f=sqrNorm(x-z)/2;
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
                    x = reshape(x,length(temp(:)),[]);
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
                    weight = xi./(sqrtSSqrMu.^3);
                    h = @(x,opt) hessian(weight,x,opt);
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
            f=sqrNorm(y-PhiAlpha)/2;
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

        function [f,g,h] = poissonModel(alpha,Phi,Phit,y)
            PhiAlpha=Phi(alpha);
            if(norm(alpha)<1e-14) PhiAlpha=ones(size(PhiAlpha)); end;
            PhiAlpha(PhiAlpha<=0)=1;
            f=sum(PhiAlpha(:))-innerProd(y,log(PhiAlpha));
            if(isnan(f)) keyboard; end
            if(nargout>=2)
                g=Phit(  1-y./PhiAlpha  );
                if(nargout>=3)
                    weight=y./(PhiAlpha.^2);
                    h=@(x,opt) hessian(weight,x,opt);
                end
            end
            function hh=hessian(weight,x,opt)
                z = Phi(x);
                if(opt==1)
                    hh = Phit(weight.*z);
                else
                    hh = innerProd(z,weight.*z);
                end
            end
        end

        % The following model assumes the y=exp(-Φα)
        function [f,g,h] = poissonModelLogLink0(alpha,Phi,Phit,y)
            PhiAlpha=Phi(alpha); weight = exp(-PhiAlpha);
            f=sum(weight)+innerProd(y,PhiAlpha);
            if(nargout>=2)
                g=Phit(  y(:) - weight  );
                if(nargout>=3)
                    h=@(x,opt) hessian(weight,x,opt);
                end
            end
            function hh=hessian(weight,x,opt)
                z = Phi(x);
                if(opt==1)
                    hh = Phit(weight.*z);
                else
                    hh = innerProd(z,weight.*z);
                end
            end
        end

        % The following model assumes the y=I_0 * exp(-Φα), where I_0 is unknown and estimated as a function of α
        function [f,g,h] = poissonModelLogLink(alpha,Phi,Phit,y)
            PhiAlpha=Phi(alpha); weight = exp(-PhiAlpha);
            sumy=sum(y); sumWeight=sum(weight);
            f=sumy*log(sumWeight)+innerProd(y,PhiAlpha);
            if(nargout>=2)
                g=Phit(  y(:) - weight*sumy/sumWeight  );
                if(nargout>=3)
                    h=@(x,opt) hessian(x,opt);
                end
            end
            function hh = hessian(x,opt)
                z = Phi(x);
                if(opt==1)
                    hh = Phit(weight.*z-weight*innerProd(weight,z)/sumWeight)*sumy/sumWeight;
                else
                    w = innerProd(weight,z);
                    hh = (innerProd(z,weight.*z)-w*w/sumWeight)*sumy/sumWeight;
                end
            end
        end

        function [dif,l1norm]=testSparsity(x)
            k=0;
            for perce=0.7:0.01:0.98
                k=k+1;
                for i=1:9
                    for j=1:4
                        level = i; wave = 2*j;
                        [dif(i,j,k),~,l1norm(i,j)]=Utils.getError(x,perce,level,wave);
                    end
                end
                minDif=dif(:,:,k); minDif=min(minDif(minDif>0));
                [i,j]=find(dif(:,:,k)==minDif);
                fprintf('level=%d,daub=%d  :  %g%% of signal achieves %g%%\n',i,2*j,(1-perce)*100,dif(i,j,k)*100);
            end
            [i,j]=find(l1norm==min(l1norm(l1norm>0)));
            fprintf('level=%d,daub=%d  :  |x|_1 of signal achieves %g\n',i,2*j,l1norm(i,j));
        end

        function [dif,coef,l1norm]=getError(x,perce,level,wave)
            h=daubcqf(wave);
            [wav,L] = mdwt(x,h,level);
            [~,idx]=sort(abs(wav(:)));
            coef=wav(idx);
            coef=coef(end:-1:1);
            %figure(2); plot(wav(idx)); hold on;
            wav(idx(1:floor(perce*length(x(:)))))=0;
            recX=midwt(wav,h,level);
            %figure(1); plot(recX); hold on; plot(x,'r');
            dif=sqrNorm(x-recX)/sqrNorm(x);
            l1norm=sum(abs(coef(:)));
        end
    end
end

