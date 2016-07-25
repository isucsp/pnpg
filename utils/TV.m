classdef TV < handle
    methods(Static,Access=private)
        function out = weight(i,j,k,l)
            %               /k+l-2\ /i+j-k-l-1\
            %               \ k-1 / \   j-l   /
            % approximate   -------------------
            %                     /i+j-2\
            %                     \ i-1 /
            % using Stirling's approximation
             
            if(i<k || j<l) 
                global strlen
                strlen=0;
                fprintf('TV.weight(i=%d, j=%d, k=%d, l=%d) needs larger i and j\n',i,j,k,l);
                return;
            end
            out = -0.5*log(2*pi);
            t=k+l-2;     out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=k-1;       out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=l-1;       out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);

            t=i+j-k-l-1; out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=i-1;       out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=j-1;       out = out + ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=i-k-1;     out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=j-l;       out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            t=i+j-2;     out = out - ((t+0.5)*log(t)+1/12/t-1/360/t^3+1/1260/t^5);
            out = exp(out); 
        end
        function out = U1(g)
            [I,J]=size(g);
            G = zeros(2*(I+1),2*(J+1));
            G(2:I+1,2:J+1)=g;

            % size of G is 2(I+1) x 2(J+1)
            % index of G is from 0 to 2(I+1)-1,   0 to 2(J+1)-1

            p1p=zeros(2*(I+1),2*(J+1));
            p1p(1:I,1)=1;
            P1=ifft2(fft2(G).*fft2(p1p));
            P1=P1(2:I+1,2:J+1);

            s=sum(P1,2);
            H = zeros(I,J);
            H(:,1:J-1)=H(:,1:J-1)-repmat(s,1,J-1);
            H(:,2:J  )=H(:,2:J  )-repmat(s,1,J-1);
            H=H-2*(I-1)*P1;
            for s=1:I
                H(s,:)=(2*s-1)*P1(I,:);
            end
            H=H/(I+J-2)/2;
            out=max(abs(H(:)));
        end
    end
    methods(Static)
        function [newX,innerSearch]=denoise(x,u,innerThresh,maxInnerItr,maskmt,tvType,lb,ub)
            if(~exist('tvType','var')) tvType='iso'; end
            if(~exist('lb','var')) lb=0; end
            if(~exist('ub','var')) ub=inf; end
            pars.print = 0;
            pars.tv = tvType;
            pars.MAXITER = maxInnerItr;
            pars.epsilon = innerThresh; 

            if(~exist('maskmt','var') || isempty(maskmt))
                mask.a=@(xx) xx(:);
                mask.b=@(xx) reshape(xx,sqrt(length(xx(:))),[]);
            else
                maskIdx=find(maskmt~=0);
                n=size(maskmt);
                mask.a=@(xx) maskFunc(xx,maskIdx);
                mask.b=@(xx) maskFunc(xx,maskIdx,n);
            end

            [newX,innerSearch]=denoise_bound_mod(mask.b(x),u,lb,ub,pars);
            newX=mask.a(newX);
        end
        function u_max = tvParUpBound(g,mask)
            global strlen
            if(~exist('mask','var') || isempty(mask))
                strlen=0;
                fprintf('\nempty mask in Utils.tvParUpBound\n');
                G=g;
            else
                maskIdx=find(mask~=0);
                n=size(mask);
                G=maskFunc(g,maskIdx,n);
            end
            [m,n]=size(G);
            u_max=+inf;
            temp=0;
            t2=0;
            for i=1:m;
                t2=t2+G(i,:);
                temp=max(temp,max(abs(t2(:))));
            end
            u_max=min(temp,u_max);
            t2=0; temp=0;
            for i=1:n;
                t2=t2+G(:,i);
                temp=max(temp,max(abs(t2(:))));
            end
            u_max=min(temp,u_max);
        end

        function o = upperBoundU(g)
            o=max(TV.U1(g),TV.U1(g'));
        end
        function o = upperBoundU_dual(g)
            [I,J]=size(g);
            h=zeros(I,J);
            h(:,1)=logspace(log10(0.5),log10(0.5)*I,I)';
            h(1,:)=logspace(log10(0.5),log10(0.5)*J,J);
            for i=2:I
                for j=2:J
                    h(i,j)=(h(i-1,j)+h(i,j-1))/2;
                end
            end
            o=conv2(h,g);
            o=abs(o(1:I,1:J));
            %figure; showImg(o); pause;
            o=max(o(:));
        end
        function o = upperBoundU_admm(g)
            [I,J]=size(g);
            p=zeros(I-1,J);
            q=zeros(I,J-1);
            lambda1=p+0/(4*I*J-2*I-2*J);
            lambda2=p+0/(4*I*J-2*I-2*J);
            lambda3=q+0/(4*I*J-2*I-2*J);
            lambda4=q+0/(4*I*J-2*I-2*J);
            %lambda1(1,1)=0.5; lambda3(1,1)=0.5;
            y=zeros(I,J);
            rho=1;
            cnt=0;

            while(cnt<100)
                cnt=cnt+1;

                u=sum([lambda1(:); lambda2(:); lambda3(:); lambda4(:)]);
                p=-inv_AtA(At(B(q)+g+y/rho)+(lambda1-lambda2)/rho);
                q=-inv_BtB(Bt(A(p)+g+y/rho)+(lambda3-lambda4)/rho);
                %keyboard;  % see value of p and q
                step=1/(4*I*J-2*I-2*J)/2;
                lambda1=max(0,lambda1+step*( p-u));
                lambda2=max(0,lambda2+step*(-p-u));
                lambda3=max(0,lambda3+step*( q-u));
                lambda4=max(0,lambda4+step*(-q-u));
                y=y+rho*(A(p)+B(q)+g);
                fprintf('lambda: %g %g step=%g, u=%g %g\n',...
                    min([lambda1(:); lambda2(:); lambda3(:); lambda4(:)]),...
                    max([lambda1(:); lambda2(:); lambda3(:); lambda4(:)]),...
                    step, max(abs([p(:); q(:)])), u);

                %if(cnt>10) % prevent excessive back and forth adjusting
                %    if(difS>10*residual)
                %        rho=rho/2 ; cnt=0;
                %    elseif(difS<residual/10)
                %        rho=rho*2 ; cnt=0;
                %    end
                %end
            end
            o=max(abs([p(:); q(:)]));

            function x = A(p)
                [I,J]=size(p); I=I+1;
                x=zeros(I,J);
                x(1:I-1,:)=p;
                x(2:I,:)=x(2:I,:)-p;
            end
            function p = At(x)
                [I,J]=size(x);
                p=x(1:I-1,:)-x(2:I,:);
            end
            function invP = inv_AtA(p)
                [I,J]=size(p); I=I+1;
                [k,i]=meshgrid(1:I-1);
                matrix=min(k,i).*min(I-k,I-i);
                invP=matrix*p/I;
            end
            function x = B(q)
                [I,J]=size(q); J=J+1;
                x=zeros(I,J);
                x(:,1:J-1)=q;
                x(:,2:J)=x(:,2:J)-q;
            end
            function q = Bt(x)
                [I,J]=size(x);
                q=x(:,1:J-1)-x(:,2:J);
            end
            function invQ = inv_BtB(q)
                [I,J]=size(q); J=J+1;
                [k,j]=meshgrid(1:J-1);
                matrix=min(k,j).*min(J-k,J-j);
                invQ=q*matrix/J;
            end
        end

        function o = upperBoundU_admm2(g)
            [I,J]=size(g);
            p=zeros(I-1,J);
            q=zeros(I,J-1);
            y=zeros(I,J);
            u=1;
            rho=1;

            cnt=0;
            ii=0;

            opts.printEvery=inf;
            opts.maxTotalIts=5e5;
            opts.maxIts=1000;
            opts.factr=1e-8/eps;
            normG=norm(g,'fro');
            relDif=1;

            while(relDif>=1e-3)
                ii=ii+1; cnt=cnt+1;

                preP=p; preQ=q; preU=u;

                % objective: 0.5*||Ap+Bq+g/u+y/rho/u||_F^2
                % subject to the [-1,1] box
                opts.x0=[p(:); q(:)];
                [x0,cost,info]=lbfgsb(@(x) quadbox(x,(rho*g+y)/(rho*u),@A,@At,@B,@Bt),...
                    -1*ones(2*I*J-I-J,1),1*ones(2*I*J-I-J,1), opts);
                p=reshape(x0(1:(I-1)*J),I-1,J);
                q=reshape(x0((I-1)*J+1:end),I,J-1);
                cost=cost;
                numItr=info.iterations;

                newG=A(p)+B(q);
                u=-(1+sum((rho*g(:)+y(:)).*newG(:)))/(1+rho*norm(newG,'fro')^2);

                y=y+rho*(u*newG+g);

                residual=norm(u*newG+g,'fro');
                difPQ=norm((preU-u)*newG,'fro');
                gap=y(:)'*reshape(u*A(p)+u*B(q)+g,[],1);

                relDif=max(residual, difPQ)/normG;

                fprintf('itr=%d, u=%g relDif=%g gap=%g\n',ii, u, relDif, gap);
                %fprintf('itr=%d, u=%g %g residual=%g rho=%g difPQ=%g, gap=%g numItr=%d normG=%g\n',...
                %   ii, max(abs([p(:); q(:)])), u, residual, rho, difPQ,gap, numItr, normG);

                %if(cnt>10) % prevent excessive back and forth adjusting
                %    if(difPQ>10*residual)
                %        rho=rho/2 ; cnt=0;
                %    elseif(difPQ<residual/10)
                %        rho=rho*2 ; cnt=0;
                %    end
                %end
            end
            o=u;

            function x = A(p)
                [I,J]=size(p); I=I+1;
                x=zeros(I,J);
                x(1:I-1,:)=p;
                x(2:I,:)=x(2:I,:)-p;
            end
            function p = At(x)
                [I,J]=size(x);
                p=x(1:I-1,:)-x(2:I,:);
            end
            function invP = inv_AtA(p)
                [I,J]=size(p); I=I+1;
                [k,i]=meshgrid(1:I-1);
                matrix=min(k,i).*min(I-k,I-i);
                invP=matrix*p/I;
            end
            function x = B(q)
                [I,J]=size(q); J=J+1;
                x=zeros(I,J);
                x(:,1:J-1)=q;
                x(:,2:J)=x(:,2:J)-q;
            end
            function q = Bt(x)
                [I,J]=size(x);
                q=x(:,1:J-1)-x(:,2:J);
            end
            function invQ = inv_BtB(q)
                [I,J]=size(q); J=J+1;
                [k,j]=meshgrid(1:J-1);
                matrix=min(k,j).*min(J-k,J-j);
                invQ=q*matrix/J;
            end

            function [f,g] = quadbox(x,c,A,At,B,Bt)
                % Err= z-f(theta)
                [I,J]=size(c);
                p=reshape(x(1:(I-1)*J),I-1,J);
                q=reshape(x((I-1)*J+1:end),I,J-1);
                r=A(p)+B(q)+c;
                f=0.5*norm(r,'fro')^2;

                if(nargout>1)
                    g=[reshape(At(r),[],1); reshape(Bt(r),[],1)];
                end
            end

        end
    end
end

