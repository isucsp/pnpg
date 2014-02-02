classdef Spline
    % Find the case when active set wanders
    properties
        polyOut
        plot
    end 

    methods
        function obj = Spline(sType,kappa)
            % Method help here
            if(nargin>0)
                if(nargin<2)
                    kappa = logspace(log10(0.03),log10(30),100);
                end
                switch lower(sType)
                    case 'dis'
                    case 'b0'
                        obj.polyOut = @(sss,III) b0Iout(kappa,sss,III);
                    case 'b1'
                end
            end
        end
    end

    methods(Static)
        function genEquations(sType)
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
                    funcpz = subs(diff(simple(funcp*s^3),s,3),s,0)/6
                    funcppz = subs(diff(simple(funcpp*s^4),s,4),s,0)/24
            end
            out.func = func; out.funcp = funcp; out.funcpp = funcpp;
            out.funcz = funcz; out.funcpz = funcpz; out.funcppz = funcppz;
        end

        function [BL,sBL,ssBL] = b0Iout(kappa, s, I)
            % kappa is a column vector, representing the nodes
            % for b0-spline, length(I)=length(kappa)-1;
            % B-0 spline with nodes be kappa

            if(nargin==0)       % start to test
                kappa = logspace(log10(0.03),log10(30),100);
                s = -0.1:0.001:10;
                [BL, sBL, ssBL] = b0Iout1(kappa,s);
                return;
            end

            kappa = kappa(:); s = s(:);
            idx = find(abs(s)<=eps);
            if(nargin>=3)
                BL = zeros(length(s),1);
            else
                BL = zeros(length(s),length(kappa)-1);
            end
            if(nargout>1) sBL = BL; end
            if(nargout>2) ssBL = BL; end
            expksR = exp(-kappa(1)*s);
            for i=1:length(kappa)-1
                expksL = expksR;
                expksR = exp(-kappa(i+1)*s);
                temp = (expksL-expksR)./s;
                temp(idx) = kappa(i+1)-kappa(i);
                if(nargin>=3) BL = BL + temp*I(i);
                else BL(:,i) = temp; end

                if(nargout>1)
                    temp = ((s*kappa(i)+1).*expksL-(s*kappa(i+1)+1).*expksR)./s./s;
                    temp(idx) = (kappa(i+1)^2-kappa(i)^2)/2;
                    if(nargin>=3) sBL = sBL + temp*I(i);
                    else sBL(:,i) = temp; end
                end
                if(nargout>2)
                    temp = ( ((kappa(i)*s+2)*kappa(i).*s+2).*expksL ...
                        -((kappa(i+1)*s+2)*kappa(i+1).*s + 2).*expksR )./s./s./s;
                    temp(idx) = (kappa(i+1)^3-kappa(i)^3)/3;
                    if(nargin>=3) ssBL = ssBL + temp*I(i);
                    else ssBL(:,i) = temp; end
                end
            end
        end
    end

end
