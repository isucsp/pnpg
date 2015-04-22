function [uk,ui] = upiota(epsilon, kappa, iota)
    if(nargin==0)
        symbol='U';
        [epsilon,iota]=readSpectrum('tungsten',140,0);
        loadXrayMassCoef
        for i=1:length(symbols)
            if(strcmpi(symbol,symbols{i}))
                massAttenCoef=mac{i};
                clear 'mac' 'comp' 'symbols' 'zaid'
                break;
            end
        end
        massAttenCoef(:,1)=massAttenCoef(:,1)*1e3;
        idx = massAttenCoef(:,1)>=min(epsilon)...
            & massAttenCoef(:,1)<=max(epsilon);
        h= max(1,find(idx, 1 )-4);
        t= min(length(idx),find(idx, 1, 'last' )+4);
        massAttenCoef=massAttenCoef(h:t,:);

        temp=massAttenCoef(2:end,1)-massAttenCoef(1:end-1,1);
        jump = find(temp==0);
        head=[1; jump(:)+1]; tail=[jump(:); length(massAttenCoef(:,1))];
        for i=1:length(head)
            h=head(i); t=tail(i);
            idx=find(epsilon>=massAttenCoef(h,1)&epsilon<=massAttenCoef(t,1));
            kappa(idx)=exp(interp1(log(massAttenCoef(h:t,1)),...
                log(massAttenCoef(h:t,2)), log(epsilon(idx)),'spline'));
        end
        kappa=kappa(:);

        [uk,ui]=upiota(epsilon, kappa, iota);

        if(nargout==0)
            figure; semilogy(epsilon,kappa,'r-'); hold on;
            loglog(massAttenCoef(:,1),massAttenCoef(:,2),'b*');
            title('\kappa(\epsilon)');
            legend('extracted \kappa(\epsilon)', 'original data');
            xlabel('\epsilon (keV)'); ylabel('\kappa (cm^2/g)');

            figure; plot(epsilon,iota,'-');
            title('spectrum tungsten 140keV with 0 relative voltage ripple');
            xlabel('\epsilon (keV)');

            figure; semilogx(uk,ui);
            title('\iota(\kappa)');
            xlabel('\kappa (cm^2/g)');
        end
        return;
    end
    uk=logspace(log10(min(kappa)),log10(max(kappa)),length(kappa))';
    temp=kappa(2:end)-kappa(1:end-1);
    jump = find(temp>0);
    head=[1; jump(:)+1]; tail=[jump(:); length(kappa)];
    ui=0*uk;
    for i=1:length(head)
        h=head(i); t=tail(i);
        temp1 = [epsilon(h); (epsilon(h:t-1)+epsilon(h+1:t))/2; epsilon(t)];
        temp2 = [  kappa(h); (  kappa(h:t-1)+  kappa(h+1:t))/2;   kappa(t)];
        temp1 = temp1(2:end)-temp1(1:end-1);
        temp2 = temp2(2:end)-temp2(1:end-1);
        idx=(uk>=kappa(t) & uk<=kappa(h));

        temp=abs(iota(h:t).*temp1./temp2);
        ui(idx)=ui(idx)+interp1(log(kappa(h:t)),temp,log(uk(idx)),'spline');
    end
end

