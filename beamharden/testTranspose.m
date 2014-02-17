function R=testTranspose(A,At,N,m,strA,NN);
    %% A is a matrix Nxm
    times=10;
    if(nargin>5) times=NN; end

    c=randperm(m); c=c(1:times);
    x1=zeros(m,1);
    y2=zeros(N,1);
    str='';
    fprintf('Verifying Symmetry of %s and %st... ',strA,strA);
    for i=1:length(c)
        x1(c(i))=1;
        y1=A(x1);
        x1(c(i))=0;
        tempIdx=find(y1~=0);
        if(isempty(tempIdx)) continue; end
        temp=randperm(length(tempIdx));
        r(i)=tempIdx(temp(1));

        y2(r(i))=1;
        x2=At(y2);
        y2(r(i))=0;

        R(i,:)=[y1(r(i)),x2(c(i))];
        temp=length(str);
        str=sprintf('%d/%d',i,times);
        fprintf([repmat('\b',1,temp) '%s'],str);
    end
    fprintf('\nerror=%g%%\n',norm(R(:,1)-R(:,2))/norm(R(:,1))*100);
