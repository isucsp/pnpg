
pp=0;
while(1)
    pp=pp+1;
    zhz=Z'*h0*Z; temp=min(eig(zhz));
    if(temp<eps)
        if(-temp>minZHZ) minZHZ=-temp; end
        zhz=zhz+minZHZ*eye(size(zhz));
    end

    deltaIe=zhz\(Z'*diff0); deltaIe=Z*deltaIe;

    estOptG=diff0-h0*deltaIe;
    temp=B(Q,:)*estOptG;
    temp=(B(Q,:)*B(Q,:)')\temp;
    lambda=ones(size(Q)); lambda(Q)=temp;
    
    %lambda=0;

    if(any(lambda<0))
        temp=find(lambda<0);
        [~,temp1]=sort(abs(temp-1/2-E/2),'descend');
        temp=temp(temp1);
        %Q(temp(end))=false;
        Q = Q & (lambda~=min(lambda));
        fprintf('%s\n', char(Q(:)'+'0') );
        Z = null(B(Q,:),'r');
        asIdx=asIdx+1;
        ASactive.itr(asIdx)=p;
        ASactive.Ie{asIdx}=Ie;
    else
        % determine the maximum possible step size
        constrainMargin = B*Ie-b;
        if(any(constrainMargin<-eps))
            display(Ie);
            %error('current Ie violate B*I>=b constraints');
        end
        temp = B*deltaIe;
        temp1 = inf*ones(size(temp));
        temp1(temp>1e-14) = constrainMargin(temp>1e-14)./temp(temp>1e-14);
        maxStep = min( temp1 );
        collide = ((temp>1e-14) & (temp1==maxStep));
        if(maxStep>eps)
            break;
        else
            Q = (Q | collide);fprintf('%s\n', char(Q(:)'+'0') );
            Z = null(B(Q,:),'r');
        end
    end
end


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

Ie(Ie<eps)=0;
if(B(end,:)*Ie<b(end)) Ie=b(end)/(B(end,:)*Ie)*Ie; end
if(stepSz==maxStep)
    Q = (Q | collide);
    fprintf( '%s\n', char(Q(:)'+'0') );
    Z = null(B(Q,:),'r');
    asIdx=asIdx+1;
    ASactive.itr(asIdx)=p;
    ASactive.Ie{asIdx}=Ie;
end

out.llI(p)=newCost; 
%if(out.stepSz<1e-10) break; end
if(deltaNormIe<thresh2) IeReady=1; break; end


