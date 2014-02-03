classdef Spline < handle
    % Find the case when active set wanders
    properties
        sType
        kappa
        polyIout
        plotSpectrum
    end 

    methods
        function obj = Spline(sType,kappa)
            % Method help here
            if(nargin>0)
                obj.sType = sType;
                if(nargin<2)
                    kappa = logspace(log10(0.03),log10(30),100);
                end
                obj.kappa = kappa;
                switch lower(sType)
                    case 'dis'
                        obj.polyIout = @(sss,III) Spline.disIout(kappa,sss,III);
                    case 'b0'
                        obj.polyIout = @(sss,III) Spline.b0Iout(kappa,sss,III);
                    case 'b1'
                        obj.polyIout = @(sss,III) Spline.b1Iout(kappa,sss,III);
                end
            end
        end
        function setPlot(obj, trueMu, trueIe, epsilon)
            switch lower(obj.sType)
                case 'dis'
                    obj.plotSpectrum = @(III) Spline.plotDisUpiota(trueMu(1:end-1),...
                        abs(trueIe(1:end-1)...
                        .*(epsilon(2:end)-epsilon(1:end-1))...
                        ./(trueMu(2:end)-trueMu(1:end-1))), ...
                        obj.kappa, III);
                case 'b0'
                    obj.plotSpectrum = @(III) Spline.plotB0Upiota(trueMu(1:end-1),...
                        abs(trueIe(1:end-1)...
                        .*(epsilon(2:end)-epsilon(1:end-1))...
                        ./(trueMu(2:end)-trueMu(1:end-1))), ...
                        obj.kappa, III);
                case 'b1'
                    obj.plotSpectrum = @(III) Spline.plotB1Upiota(trueMu(1:end-1),...
                        abs(trueIe(1:end-1)...
                        .*(epsilon(2:end)-epsilon(1:end-1))...
                        ./(trueMu(2:end)-trueMu(1:end-1))), ...
                        obj.kappa, III);
            end
        end
    end

    methods(Static)
        function out = genEquations(sType)
            syms k s a b c
            switch lower(sType)
                case 'dis'
                case 'b0'
                    func = int(exp(-k*s),k,a,b)
                    funcp = simple(-diff(func,s))
                    funcpp = simple(diff(func,s,2))

                    funcz = subs(diff(func*s,s),s,0)
                    funcpz = subs(diff(simple(funcp*s^2),s,2),s,0)/2
                    funcppz = subs(diff(simple(funcpp*s^3),s,3),s,0)/6
                case 'b1'
                    func = int((k-a)*exp(-k*s),k,a,b)/(b-a);
                    func = simple(func + int((c-k)*exp(-k*s),k,b,c)/(c-b))
                    funcp = simple(-diff(func,s))
                    funcpp = simple(diff(func,s,2))

                    funcz = subs(diff(func*s^2,s,2),s,0)/2
                    funcz2 = subs(diff(func*s^2,s,3),s,0)/6

                    funcpz = subs(diff(simple(funcp*s^3),s,3),s,0)/6
                    funcpz2 = subs(diff(simple(funcp*s^3),s,4),s,0)/24
                    funcpz3 = subs(diff(simple(funcp*s^3),s,5),s,0)/120

                    funcppz = subs(diff(simple(funcpp*s^4),s,4),s,0)/24
                    funcppz2 = subs(diff(simple(funcpp*s^4),s,5),s,0)/120
                    funcppz3 = subs(diff(simple(funcpp*s^4),s,6),s,0)/720
                    funcppz4 = subs(diff(simple(funcpp*s^4),s,7),s,0)/720/7
            end
            out.func = func; out.funcp = funcp; out.funcpp = funcpp;
            out.funcz = funcz; out.funcpz = funcpz; out.funcppz = funcppz;
        end

        function [BL,sBL,ssBL] = b1Iout(k, s, I)
            % k is a column vector, representing the nodes
            % for b1-spline, length(I)=length(k)-1;
            % B-1 spline with nodes be k

            if(nargin==0)       % start to test
                k = logspace(log10(0.001),log10(30),100);
                s = -0.1:0.001:10;
                [BL, sBL, ssBL] = Spline.b1Iout(k,s);
                figure; semilogy(s,BL);
                figure; semilogy(s,sBL);
                figure; semilogy(s,ssBL);
                return;
            end

            k = k(:); s = s(:);
            idx = abs(s)<=eps;
            if(~isempty(I))
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(k)-1);
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end

            expksC = exp(-k(1)*s);
            expksR = exp(-k(2)*s);
            for i=1:length(k)-2
                expksL = expksC; expksC = expksR;
                expksR = exp(-k(i+2)*s);

                temp=-((k(i)-k(i+1))*expksR-(k(i)-k(i+2))*expksC-(k(i+2)-k(i+1))*expksL)...
                    ./(s.^2*(k(i)-k(i+1))*(k(i+1)-k(i+2)));

                temp(idx)=-(k(i+1)*(k(i)^2-k(i+2)^2)-k(i)*k(i+1)^2+k(i)*k(i+2)^2-k(i)^2*k(i+2)+k(i+1)^2*k(i+2))/(2*(k(i)-k(i+1))*(k(i+1)-k(i+2)));
                if(~isempty(I)) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp = -(k(i+1)*(k(i)*expksL-k(i+2)*expksR)-k(i)*k(i+1)*expksC-k(i)*k(i+2)*expksL+k(i)*k(i+2)*expksR+k(i+1)*k(i+2)*expksC)...  
                        ./(s.^2*(k(i)-k(i+1))*(k(i+1)-k(i+2)))...
                        -(2*(k(i)*expksR-k(i)*expksC-k(i+2)*expksL+k(i+2)*expksC+k(i+1)*(expksL-expksR)))...
                        ./(s.^3*(k(i)-k(i+1))*(k(i+1)-k(i+2)));

                    temp(idx)=-(k(i+1)*(k(i)^3-k(i+2)^3)-k(i)*k(i+1)^3+k(i)*k(i+2)^3-k(i)^3*k(i+2)+k(i+1)^3*k(i+2))/(6*(k(i)-k(i+1))*(k(i+1)-k(i+2)));
                    if(~isempty(I)) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp=-(4*(k(i+1)*(k(i)*expksL-k(i+2)*expksR)-k(i)*k(i+1)*expksC-k(i)*k(i+2)*expksL+k(i)*k(i+2)*expksR+k(i+1)*k(i+2)*expksC))...
                        ./(s.^3*(k(i)-k(i+1))*(k(i+1)-k(i+2)))...
                        -(k(i+1)*(k(i)^2*expksL-k(i+2)^2*expksR)-k(i)*k(i+1)^2*expksC-k(i)^2*k(i+2)*expksL+k(i)*k(i+2)^2*expksR+k(i+1)^2*k(i+2)*expksC)...
                        ./(s.^2*(k(i)-k(i+1))*(k(i+1)-k(i+2)))...
                        -(6*(k(i)*expksR-k(i)*expksC-k(i+2)*expksL+k(i+2)*expksC+k(i+1)*(expksL-expksR)))...
                        ./(s.^4*(k(i)-k(i+1))*(k(i+1)-k(i+2)));

                    temp(idx)=-(k(i+1)*(k(i)^4-k(i+2)^4)-k(i)*k(i+1)^4+k(i)*k(i+2)^4-k(i)^4*k(i+2)+k(i+1)^4*k(i+2))/(12*(k(i)-k(i+1))*(k(i+1)-k(i+2)));
                    if(~isempty(I)) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end

        function [BL,sBL,ssBL] = b0Iout(k, s, I)
            % k is a column vector, representing the nodes
            % for b0-spline, length(I)=length(k)-1;
            % B-0 spline with nodes be k

            if(nargin==0)       % start to test
                k = logspace(log10(0.001),log10(30),100);
                s = -0.1:0.001:10;
                [BL, sBL, ssBL] = Spline.b0Iout(k,s);
                figure; semilogy(s,BL);
                figure; semilogy(s,sBL);
                figure; semilogy(s,ssBL);
                return;
            end

            k = k(:); s = s(:);
            idx = abs(s)<=eps;
            if(~isempty(I))
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(k)-1);
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end
            expksR = exp(-k(1)*s);
            for i=1:length(k)-1
                expksL = expksR;
                expksR = exp(-k(i+1)*s);
                temp = (expksL-expksR)./s;
                temp(idx) = k(i+1)-k(i);
                if(~isempty(I)) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp = ((s*k(i)+1).*expksL-(s*k(i+1)+1).*expksR)./s./s;
                    temp(idx) = (k(i+1)^2-k(i)^2)/2;
                    if(~isempty(I)) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp = ( ((k(i)*s+2)*k(i).*s+2).*expksL ...
                        -((k(i+1)*s+2)*k(i+1).*s + 2).*expksR )./s./s./s;
                    temp(idx) = (k(i+1)^3-k(i)^3)/3;
                    if(~isempty(I)) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end

        function [BL,sBL,ssBL] = disIout(kappa, s, I)

            if(nargin==0)       % start to test
                kappa = logspace(log10(0.001),log10(30),100);
                s = -0.1:0.001:10;
                [BL, sBL, ssBL] = Spline.disIout(kappa,s);
                figure; semilogy(s,BL);
                figure; semilogy(s,sBL);
                figure; semilogy(s,ssBL);
                return;
            end

            kappa = kappa(:); s = s(:);
            if(~isempty(I))
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(kappa)-1);
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end

            for i=1:length(kappa)
                temp = exp(-kappa(i)*s);
                if(~isempty(I)) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp = temp*kappa(i);
                    if(~isempty(I)) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp = temp*kappa(i);
                    if(~isempty(I)) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end


        function plotB1Upiota(trueMu, trueUpiota, mu, Ie)
            loglog(trueMu,trueUpiota,'r.-'); hold on;
            loglog(mu,[0; Ie; 0],'*-'); hold off;
            %ylim([1e-10 1]);
            xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
        end

        function plotB0Upiota(trueMu, trueUpiota, mu, Ie)
            loglog(trueMu,trueUpiota,'r.-'); hold on;
            mu=reshape([mu(:)';mu(:)'],[],1);
            mu(1)=[]; mu(end)=[];
            Ie=reshape([Ie(:)';Ie(:)'],[],1);
            loglog(mu,Ie,'*-'); hold off;
            %ylim([1e-10 1]);
            xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
        end

        function plotDisUpiota(trueMu, trueUpiota, mu, Ie)
            loglog(trueMu,trueUpiota,'r.-'); hold on;
            temp=[mu(:); mu(end)^2/mu(end-1)];
            mu=reshape([temp';temp'],[],1);
            mu(1)=[]; mu(end)=[];
            temp=temp(2:end)-temp(1:end-1);
            Ie = Ie./temp;
            Ie=reshape([Ie(:)';Ie(:)'],[],1);
            loglog(mu,Ie,'*-'); hold off;
            %ylim([1e-10 1]);
            xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
        end
    end

end
