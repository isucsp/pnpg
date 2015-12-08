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
    end
end

