function R=testTranspose(A,At,N,m,strA,NN);
%testTranspose  test the symmetry of two Matrix A and At
%   Syntax:
%   R = testTranspose(A,At,N,m,strA,NN); 
%   A         function handle of A, whose size is Nxm
%   At        function handle of A'
%   strA      a string, Information to show
%   NN        number of trials to test.
%

    if(nargin==0)
        help testTranspose
    end
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
        if(abs((R(i,1)-R(i,2))/(abs(R(i,1))+abs(R(i,2))))>0.1)
            fprintf('\nBig error at i=%d, @(%d,%d) error=%s\n',i,r(i),c(i),sprintf('%g-%g=%g',R(i,1),R(i,2),R(i,1)-R(i,2)));
            str='';
        end
    end
    fprintf('\nerror=%g%%\n',norm(R(:,1)-R(:,2))/norm(R(:,1))*100);
