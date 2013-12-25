
loadMac

if(0) %verify the composite's mass attenuation coefficient
    epsilon=logspace(-3,log10(20),1024);
    epsilon=epsilon(:);
    mu=zeros(size(epsilon));
    for i=1:size(weight,1)
        temp=mac{weight(i,1)};
        delta=temp(2:end,1)-temp(1:end-1,1);
        idx=find(delta==0);
        temp(idx,1)=temp(idx,1)-eps;
        newMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
        %loglog(epsilon,newMu,'k-');
        loglog(mac{weight(i,1)}(:,1),mac{weight(i,1)}(:,2),'*');
        hold on;
        loglog(epsilon,newMu);
        hold off;
        pause();
        mu=newMu*weight(i,2)+mu;
    end
    figure;
    loglog(epsilon,mu,'k');
    hold on;
    loglog(bone(:,1),bone(:,2),'r.');
    return;
end

if(0)
    figure(3);
    baseMaterial=26;
    temp=mac{baseMaterial};
    epsilon=logspace(-3,log10(20),4096);
    epsilon=epsilon(:);
    delta=temp(2:end,1)-temp(1:end-1,1);
    idx=find(delta==0);
    temp(idx,1)=temp(idx,1)-eps;
    baseMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
    original = baseMu;
    [~, idx]=sort(baseMu,'descend');
    %epsilon=epsilon(idx);
    loglog(mac{baseMaterial}(:,1),mac{baseMaterial}(:,2),'*');
    hold on;
    loglog(epsilon,baseMu);
    hold off;

    for i=1:1 %length(density)
        temp=bone.mac; %mac{i};
        delta=temp(2:end,1)-temp(1:end-1,1);
        idx1=find(delta==0);
        temp(idx1,1)=temp(idx1,1)*0.99;
        newMu=exp(interp1(log(temp(:,1)),log(temp(:,2)),log(epsilon),'linear'));
        newMu=newMu./baseMu;
        newMu=newMu; %-min(newMu(:))+eps;

        %loglog(epsilon,newMu,'k-');
        figure(21); loglog(baseMu(idx),newMu(idx),'-'); hold on;
        figure(22); plot(baseMu(idx),newMu(idx),'-'); hold on;
        figure(23); semilogy(baseMu(idx),newMu(idx),'-'); hold on;
        figure(25); semilogx(baseMu(idx),newMu(idx),'-'); hold on;
        figure(24); semilogy(epsilon,newMu,'-'); hold on;
    end
    return;
end

if(1)
    N=length(density);
    for i=1:length(density)
        color=i/N*3;
        if(color<1)
            r=mod(color,1);
            g=0; b=0;
        elseif(color<2)
            r=1;
            g=mod(color-1,1);
            k=0;
        else
            r=1;g=1;
            b=mod(color-2,1);
        end
        figure(31); loglog(mac{i}(:,1),mac{i}(:,2),'color',[r,g,b]);
        hold on;
        figure(32); 
        loglog(mac{i}(:,1),mac{i}(:,2)*density(i),'color',[r,g,b]);
        hold on;
        figure(33);
        semilogy(mac{i}(:,1)*1000,mac{i}(:,2),'color',[r,g,b]);
        hold on;
    end
end


