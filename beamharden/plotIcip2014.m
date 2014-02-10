
style = {'b-','r-','b:','r:','b--','r--'};
for i = 1:1
    for j=1:6
        for k=1:6
            figure(1);
            eval(sprintf('semilogy(1:out2%1$d{i,j}.p,out2%1$d{i,j}.RMSE,style{k})',k));
            hold on;
            figure(2);
            eval(sprintf('semilogy(1:out2%1$d{i,j}.p,out2%1$d{i,j}.cost,style{k})',k));
            hold on;
            figure(3);
            eval(sprintf('semilogy(1:out2%1$d{i,j}.p,out2%1$d{i,j}.llAlpha,style{k})',k));
            hold on;
            figure(4);
            eval(sprintf('semilogy(1:out2%1$d{i,j}.p,out2%1$d{i,j}.penAlpha,style{k})',k));
            hold on;
        end
        pause;
        close all;
    end
end

%%

style = {'b-','r-','b:','r:','b--','r--'};
for k=1:6
    eval(sprintf('out=out2%d;',k));
    for i=1:size(out,1)
        for j=1:size(out,2)
            if(~isempty(out{i,j}))
                rse(i,j)=min(out{i,j}.RMSE);
            end
        end
    end
    eval(sprintf('rse2%d=rse;',k));
end
