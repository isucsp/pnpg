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
        function setPlot(obj, trueKappa, trueIota, epsilon)
            temp1 = [epsilon(1);(epsilon(1:end-1)+epsilon(2:end))/2;epsilon(end)];
            temp2 = [trueKappa(1);(trueKappa(1:end-1)+trueKappa(2:end))/2;trueKappa(end)];
            temp1 = temp1(2:end)-temp1(1:end-1);
            temp2 = temp2(2:end)-temp2(1:end-1);
            trueUpiota=abs(trueIota.*temp1./temp2);
            switch lower(obj.sType)
                case 'dis'
                    obj.plotSpectrum = @(III) Spline.plotDisUpiota(trueKappa,trueUpiota,...
                        obj.kappa, III);
                case 'b0'
                    obj.plotSpectrum = @(III) Spline.plotB0Upiota(trueKappa,trueUpiota,...
                        obj.kappa, III);
                case 'b1'
                    obj.plotSpectrum = @(III) Spline.plotB1Upiota(trueKappa,trueUpiota,...
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

                    funcz  = simple(subs(diff(simple(func  *s  ),s,1),s,0))
                    funcz2 = simple(subs(diff(simple(func  *s  ),s,2),s,0)/2)
                    funcz3 = simple(subs(diff(simple(func  *s  ),s,3),s,0)/6)
                    funcz4 = simple(subs(diff(simple(func  *s  ),s,4),s,0)/24)
                    
                    funcpz = simple(subs(diff(simple(funcp *s^2),s,2),s,0)/2)
                    funcpz2= simple(subs(diff(simple(funcp *s^2),s,3),s,0)/6)
                    funcpz3= simple(subs(diff(simple(funcp *s^2),s,4),s,0)/24)
                    funcpz4= simple(subs(diff(simple(funcp *s^2),s,5),s,0)/120)
                    
                    funcppz= simple(subs(diff(simple(funcpp*s^3),s,3),s,0)/6)
                    funcppz2=simple(subs(diff(simple(funcpp*s^3),s,4),s,0)/24)
                    funcppz3=simple(subs(diff(simple(funcpp*s^3),s,5),s,0)/120)
                    funcppz4=simple(subs(diff(simple(funcpp*s^3),s,6),s,0)/720)
                case 'b1'
                    func = int((k-a)*exp(-k*s),k,a,b)/(b-a);
                    func =-(exp(-a*s) - exp(-b*s))/(s^2*(a - b))
                    func1 = (exp(-b*s) - exp(-c*s))/(s^2*(b - c));
                    test=simple(func+func1-int((k-a)*exp(-k*s),k,a,b)/(b-a)-int((c-k)*exp(-k*s),k,b,c)/(c-b));
                    funcp = simple(-diff(func,s))
                    funcpp = simple(diff(func,s,2))
                    latex(func);
                    latex(funcp);
                    latex(funcpp);

                    funcz  = simple(subs(diff(simple(func  *s^2),s,2),s,0)/2)
                    funcz2 = simple(subs(diff(simple(func  *s^2),s,3),s,0)/6)
                    funcz3 = simple(subs(diff(simple(func  *s^2),s,4),s,0)/24)
                    funcz4 = simple(subs(diff(simple(func  *s^2),s,5),s,0)/120)

                    funcpz = simple(subs(diff(simple(funcp *s^3),s,3),s,0)/6)
                    funcpz2= simple(subs(diff(simple(funcp *s^3),s,4),s,0)/24)
                    funcpz3= simple(subs(diff(simple(funcp *s^3),s,5),s,0)/120)
                    funcpz4= simple(subs(diff(simple(funcp *s^3),s,6),s,0)/720)

                    funcppz= simple(subs(diff(simple(funcpp*s^4),s,4),s,0)/24)
                    funcppz2=simple(subs(diff(simple(funcpp*s^4),s,5),s,0)/120)
                    funcppz3=simple(subs(diff(simple(funcpp*s^4),s,6),s,0)/720)
                    funcppz4=simple(subs(diff(simple(funcpp*s^4),s,7),s,0)/720/7)

                    func = simple(int((c-k)*exp(-k*s),k,b,c)/(c-b));
                    func = (exp(-b*s) - exp(-c*s))/(s^2*(b - c))
                    funcp = simple(-diff(func,s))
                    funcpp = simple(diff(func,s,2))
                    latex(func);
                    latex(funcp);
                    latex(funcpp);

                    funcz  = simple(subs(diff(simple(func  *s^2),s,2),s,0)/2)
                    funcz2 = simple(subs(diff(simple(func  *s^2),s,3),s,0)/6)
                    funcz3 = simple(subs(diff(simple(func  *s^2),s,4),s,0)/24)
                    funcz4 = simple(subs(diff(simple(func  *s^2),s,5),s,0)/120)

                    funcpz = simple(subs(diff(simple(funcp *s^3),s,3),s,0)/6)
                    funcpz2= simple(subs(diff(simple(funcp *s^3),s,4),s,0)/24)
                    funcpz3= simple(subs(diff(simple(funcp *s^3),s,5),s,0)/120)
                    funcpz4= simple(subs(diff(simple(funcp *s^3),s,6),s,0)/720)

                    funcppz= simple(subs(diff(simple(funcpp*s^4),s,4),s,0)/24)
                    funcppz2=simple(subs(diff(simple(funcpp*s^4),s,5),s,0)/120)
                    funcppz3=simple(subs(diff(simple(funcpp*s^4),s,6),s,0)/720)
                    funcppz4=simple(subs(diff(simple(funcpp*s^4),s,7),s,0)/720/7)

            end
            out.func = func; out.funcp = funcp; out.funcpp = funcpp;
            out.funcz = funcz; out.funcpz = funcpz; out.funcppz = funcppz;
        end

        function [BL,sBL,ssBL] = b1Iout(k, s, I)
            % k is a column vector, representing the nodes
            % for b1-spline, length(I)=length(k)-1;
            % B-1 spline with nodes be k

            % The accuracy will be at least 100*EPS% of the value at s=0
            EPS=1e-6;
            if(nargin==0)       % start to test
                k = logspace(log10(0.001),log10(30),100);
                s = linspace(-1e-2,1e-2,1000);
                [BL, sBL, ssBL] = Spline.b1Iout(k,s,[]);
                figure; semilogy(s,BL);
                figure; semilogy(s,sBL);
                figure; semilogy(s,ssBL);
                return;
            end

            k = k(:); s = s(:);
            if(~isempty(I))
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(k)-2);
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end

            expksC = exp(-k(1)*s);
            expksR = exp(-k(2)*s);
            s2=s.^2; s3=s.^3; s4=s.^4; abss=abs(s);
            for i=1:length(k)-2
                a=k(i); b=k(i+1); c=k(i+2);
                expksL = expksC; expksC = expksR; expksR = exp(-c*s);

                temp=-(expksL-expksC)./(s2*(a-b))...
                    + (expksC-expksR)./(s2*(b-c));
                temp1=[
                -a/2+c/2;
                (a^3-b^3)/(6*(a-b))-(b^3-c^3)/(6*(b-c));
                -(a^4-b^4)/(24*(a-b))+(b^4-c^4)/(24*(b-c));
                (a^5-b^5)/(120*(a-b))-(b^5-c^5)/(120*(b-c));
                ];
                % The accuracy will be at least 100*EPS% of the value at s=0
                temp2=abs(EPS*temp1(1)/temp1(end))^(1/(length(temp1)-1));
                idx = abss<temp2;
                temp(idx) = temp1(1)+temp1(2)*s(idx);

                if(~isempty(I)) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp=-(2*(expksL-expksC))./(s3*(a-b))-(a*expksL-b*expksC)./(s2*(a-b))...
                        +(2*(expksC-expksR))./(s3*(b-c))+(b*expksC-c*expksR)./(s2*(b-c));
                    temp1=[
                    -(a^3-b^3)/(6*(a-b))+(b^3-c^3)/(6*(b-c));
                    (a^4-b^4)/(12*(a-b))-(b^4-c^4)/(12*(b-c));
                    -(a^5-b^5)/(40*(a-b))+(b^5-c^5)/(40*(b-c));
                    (a^6-b^6)/(180*(a-b))-(b^6-c^6)/(180*(b-c));
                    ];
                    temp2=abs(EPS*temp1(1)/temp1(end))^(1/(length(temp1)-1));
                    idx = abss<temp2;
                    temp(idx) = temp1(1)+temp1(2)*s(idx)+temp1(3)*s2(idx);
                    if(~isempty(I)) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp=-(6*expksL-6*expksC+s2.*(a^2*expksL-b^2*expksC)+s.*(4*a*expksL-4*b*expksC))./(s4*(a-b))...
                        +(6*expksC-6*expksR+s2.*(b^2*expksC-c^2*expksR)+s.*(4*b*expksC-4*c*expksR))./(s4*(b-c));
                    
                    temp1=[
                    -((a+b)*(a^2+b^2))/12+((b+c)*(b^2+c^2))/12;
                    (a^4+a^3*b+a^2*b^2+a*b^3+b^4)*20^(-1)+(-b^4-b^3*c-b^2*c^2-b*c^3-c^4)*20^(-1);
                    -((a+b)*(a^2+a*b+b^2)*(a^2-a*b+b^2))/60+((b+c)*(b^2+b*c+c^2)*(b^2-b*c+c^2))/60;
                    (a^6+a^5*b+a^4*b^2+a^3*b^3+a^2*b^4+a*b^5+b^6)*252^(-1)+(-b^6-b^5*c-b^4*c^2-b^3*c^3-b^2*c^4-b*c^5-c^6)*252^(-1)
                    ];
                    temp2=abs(EPS*temp1(1)/temp1(end))^(1/(length(temp1)-1));
                    idx = abss<temp2;
                    temp(idx) = temp1(1)+temp1(2)*s(idx)...
                        +temp1(3)*s2(idx)+temp1(4)*s3(idx);
                    
                    if(~isempty(I)) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end

        function [BL,sBL,ssBL] = b0Iout(k, s, I)
            % k is a column vector, representing the nodes
            % for b0-spline, length(I)=length(k)-1;
            % B-0 spline with nodes be k

            EPS = 1e-6;
            if(nargin==0)       % start to test
                k = logspace(log10(0.001),log10(30),100);
                s = linspace(-1e-2,1e-2,1000);
                [BL, sBL, ssBL] = Spline.b0Iout(k,s,[]);
                figure; semilogy(s,BL);
                figure; semilogy(s,sBL);
                figure; semilogy(s,ssBL);
                return;
            end

            k = k(:); s = s(:);
            if(~isempty(I))
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(k)-1);
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end
            expksR = exp(-k(1)*s);
            s2=s.^2; s3=s.^3; abss=abs(s);
            for i=1:length(k)-1
                a=k(i); b=k(i+1);
                expksL = expksR; expksR = exp(-b*s);
                temp = (expksL-expksR)./s;
                temp1=[b-a;a^2/2-b^2/2;b^3/6-a^3/6;a^4/24-b^4/24;];
                temp2=abs(EPS*temp1(1)/temp1(end))^(1/(length(temp1)-1));
                idx = abss<temp2;
                temp(idx) = temp1(1)+temp1(2)*s(idx);
                if(~isempty(I)) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp = ((s*a+1).*expksL-(s*b+1).*expksR)./s2;
                    temp1=[b^2/2-a^2/2;a^3/3-b^3/3;b^4/8-a^4/8;a^5/30-b^5/30;];
                    temp2=abs(EPS*temp1(1)/temp1(end))^(1/(length(temp1)-1));
                    idx = abss<temp2;
                    temp(idx) = temp1(1)+temp1(2)*s(idx)+temp1(3)*s2(idx);
                    if(~isempty(I)) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp = ( ((a*s+2)*a.*s+2).*expksL ...
                        -((b*s+2)*b.*s + 2).*expksR )./s3;
                    temp1=[(b^3-a^3)/3;(a^4-b^4)/4;(b^5-a^5)/10;(a^6-b^6)/36];
                    temp2=abs(EPS*temp1(1)/temp1(end))^(1/(length(temp1)-1));
                    idx = abss<temp2;
                    temp(idx) = temp1(1)+temp1(2)*s(idx)...
                        +temp1(3)*s2(idx)+temp1(4)*s3(idx);
                    if(~isempty(I)) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end

        function [BL,sBL,ssBL] = disIout(kappa, s, I)

            if(nargin==0)       % start to test
                kappa = logspace(log10(0.001),log10(30),100);
                s = -0.1:0.001:10;
                [BL, sBL, ssBL] = Spline.disIout(kappa,s,[]);
                figure; semilogy(s,BL);
                figure; semilogy(s,sBL);
                figure; semilogy(s,ssBL);
                return;
            end

            dKappa=(kappa(3:end)-kappa(1:end-2))/2;
            kappa = kappa(2:end-1); kappa = kappa(:); s = s(:);
            if(~isempty(I))
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(kappa));
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end

            for i=1:length(kappa)
                temp = exp(-kappa(i)*s)*dKappa(i);
                if(~isempty(I)) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp = temp*(-kappa(i));
                    if(~isempty(I)) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp = temp*(-kappa(i));
                    if(~isempty(I)) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end

        function [trueMu,trueUpiota,mu,Ie]= plotB1Upiota(trueMu, trueUpiota, mu, Ie)
            Ie=[0; Ie; 0];
            if(nargout==0)
                loglog(trueMu,trueUpiota,'r.-'); hold on;
                loglog(mu,Ie,'*-'); hold off;
                %ylim([1e-10 1]);
                xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
            end
        end

        function [trueMu,trueUpiota,mu,Ie]=plotB0Upiota(trueMu, trueUpiota, mu, Ie)
            mu=reshape([mu(:)';mu(:)'],[],1);
            mu(1)=[]; mu(end)=[];
            Ie=reshape([Ie(:)';Ie(:)'],[],1);
            if(nargout==0)
                loglog(trueMu,trueUpiota,'r.-'); hold on;
                loglog(mu,Ie,'*-'); hold off;
                %ylim([1e-10 1]);
                xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
            end
        end

        function [trueMu,trueUpiota,mu,Ie]=plotDisUpiota(trueMu,trueUpiota,mu,Ie)
            mu=(mu(1:end-1)+mu(2:end))/2;
            mu=reshape([mu(:)';mu(:)'],[],1);
            mu(1)=[]; mu(end)=[];
            Ie=reshape([Ie(:)';Ie(:)'],[],1);
            if(nargout==0)
                loglog(trueMu,trueUpiota,'r.-'); hold on;
                loglog(mu,Ie,'*-'); hold off;
                %ylim([1e-10 1]);
                xlim([min(min(trueMu),mu(1)) max(max(trueMu),mu(end))]);
            end
        end
    end
end

