
while(1)
    zhz=Z'*h0*Z; temp=min(eig(zhz));
    if(temp<0)
        if(-temp>minZHZ) minZHZ=-temp; end
        zhz=zhz+minZHZ*eye(size(zhz));
    end

    deltaIe=zhz\(Z'*diff0); deltaIe=Z*deltaIe;

    estOptG=diff0-h0*deltaIe;
    temp=AA(G,:); temp=(temp*temp')\temp;
    lambda=zeros(E+1,1); lambda(G)=temp*estOptG;
    %lambda1=(AA(G,:)*inv(h0)*AA(G,:)')\(AA(G,:)*(inv(h0)*diff0));
    %lambda1'-lambda(G)'
    %deltaIe1=inv(h0)*(diff0-AA(G,:)'*lambda1);

    if(any(lambda<0))
        %temp=find(lambda==min(lambda));
        temp=find(lambda<0);
        [~,temp2]=sort(abs(temp-1/2-E/2),'descend');
        temp=temp(temp2);
        temp1=zeros(E+1,1); temp1(temp(end))=1;
        AAIe=AAIe&(temp1==0);
        [P,Pbar,Z,G,Gbar]=updateASIe(AAIe);
        asIdx=asIdx+1;
        ASactive.itr(asIdx)=p;
        ASactive.Ie{asIdx}=Ie;
    else break; end
end

temp=find(deltaIe>0);
if(isempty(temp)) opt.maxStep1=1.1;
else step=Ie(temp)./deltaIe(temp);
    temp1=find(step==min(step));
    if(min(step)==-100)
        temp1=find( step==0 & ...  
            deltaIe(temp)==max(deltaIe(temp(temp1))) );
    end
    temp=temp(temp1); [~,temp1]=sort(abs(temp-E/2-1/2),'descend');
    temp=temp(temp1(1));
    temp1=zeros(E+1,1); temp1(temp)=1;
    opt.maxStep1=Ie(temp)/deltaIe(temp);
end
opt.maxStep1=max(opt.maxStep1,0);

if(sum(deltaIe)<0)
    if(1-sum(Ie)>0)
        opt.maxStep2=(sum(Ie)-1)/sum(deltaIe);
    else opt.maxStep2=0;
    end
else opt.maxStep2=1;
end
opt.maxStep2=max(opt.maxStep2,0);
opt.maxStep =min(opt.maxStep1,opt.maxStep2);

opt.g=diff0;
dummyf=@(x) 0;
out = lineSearch(Ie,deltaIe,f0,dummyf,1,costA,opt);
deltaNormIe=out.delta; stepSzIe=out.stepSz;
if(stepSzIe==opt.maxStep)
    if(stepSzIe==opt.maxStep2)
        AAIe(end)=true;
    else
        AAIe=AAIe | (temp1==1);
    end
    [P,Pbar,Z,G,Gbar]=updateASIe(AAIe);
    asIdx=asIdx+1;
    ASactive.itr(asIdx)=p;
    ASactive.Ie{asIdx}=out.x;
end

res2(p)=out.costA; Ie=out.x; zmf=out.zmf; Ie(Ie<0)=0;
%if(out.stepSz<1e-10) break; end
if(deltaNormIe<thresh2) IeReady=1; break; end
