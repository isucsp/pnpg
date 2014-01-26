
while(1)
    zhz=Z'*h0*Z; temp=min(eig(zhz));
    if(temp<0)
        if(-temp>minZHZ) minZHZ=-temp; end
        zhz=zhz+minZHZ*eye(size(zhz));
    end

    deltaIe=zhz\(Z'*diff0); deltaIe=Z*deltaIe;

    estOptG=diff0-h0*deltaIe;
    temp=B(Q,:)*estOptG;
    temp=(B(Q,:)*B(Q,:)')\temp;
    lambda=zeros(E+1,1); lambda(Q)=temp;

    if(any(lambda<0))
        %temp=find(lambda==min(lambda));
        temp=find(lambda<0);
        [~,temp2]=sort(abs(temp-1/2-E/2),'descend');
        temp=temp(temp2);
        temp1=zeros(E+1,1); temp1(temp(end))=1;
        Q=Q&(temp1==0);
            Z = null(B(Q,:));
        asIdx=asIdx+1;
        ASactive.itr(asIdx)=p;
        ASactive.Ie{asIdx}=Ie;
    else break; end
end

% determine the maximum possible step size
constrainMargin = B*Ie-b;
if(any(constrainMargin<0))
    display(Ie);
    error('current Ie violate B*I>=b constraints');
end
temp = B*deltaIe;
maxStep = min( constrainMargin(temp>0)./temp(temp>0) );

% begin line search
pp=0; stepSz=min(1,maxStep);
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
    Q = (B*Ie-b<1e-14);
    Z = null(B(Q,:));
    asIdx=asIdx+1;
    ASactive.itr(asIdx)=p;
    ASactive.Ie{asIdx}=Ie;
end

Ie(Ie<0)=0;
out.llI(p)=newCost; 
%if(out.stepSz<1e-10) break; end
if(deltaNormIe<thresh2) IeReady=1; break; end


