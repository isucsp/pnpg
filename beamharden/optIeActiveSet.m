
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
if(isempty(temp)) maxStep1=1.1;
else step=Ie(temp)./deltaIe(temp);
    temp1=find(step==min(step));
    if(min(step)==-100)
        temp1=find( step==0 & ...  
            deltaIe(temp)==max(deltaIe(temp(temp1))) );
    end
    temp=temp(temp1); [~,temp1]=sort(abs(temp-E/2-1/2),'descend');
    temp=temp(temp1(1));
    temp1=zeros(E+1,1); temp1(temp)=1;
    maxStep1=Ie(temp)/deltaIe(temp);
end
maxStep1=max(maxStep1,0);

if(sum(deltaIe)<0)
    if(1-sum(Ie)>0)
        maxStep2=(sum(Ie)-1)/sum(deltaIe);
    else maxStep2=0;
    end
else maxStep2=1;
end
maxStep2=max(maxStep2,0);
maxStep =min(maxStep1,maxStep2);

% begin line search
pp=0; maxStep=maxStep;
stepSz=min(1,maxStep);
while(1)
    pp=pp+1;
    newX=Ie-stepSz*deltaIe;
    %newX(newX<0)=0; % force it be positive;

    [newCost,zmf]=llI(A,newX);

    if(newCost <= costA - stepSz/2*diff0'*deltaIe)
        break;
    else
        if(pp>10) stepSz=0; else stepSz=stepSz*stepShrnk; end
    end
end
% end of line search
deltaNormIe=diff0'*deltaIe;
Ie = newX;

if(stepSz==maxStep)
    if(stepSz==maxStep2)
        AAIe(end)=true;
    else
        AAIe=AAIe | (temp1==1);
    end
    [P,Pbar,Z,G,Gbar]=updateASIe(AAIe);
    asIdx=asIdx+1;
    ASactive.itr(asIdx)=p;
    ASactive.Ie{asIdx}=Ie;
end

Ie(Ie<0)=0;
out.llI(p)=newCost; 
%if(out.stepSz<1e-10) break; end
if(deltaNormIe<thresh2) IeReady=1; break; end


