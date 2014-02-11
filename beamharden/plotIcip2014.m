
style1 = {'b-','r-','b:','r:','b--','r--'};
style2 = {'b-*','r-*','b:o','r:o','b--s','r--s'};

if(0)
    for k=1:6
        eval(sprintf('out=out1%d;',k));
        rse = zeros(size(out));
        for i=1:size(out,1)
            for j=1:size(out,2)
                if(~isempty(out{i,j}))
                    figure((i-1)*size(out,2)+j);
                    loglog(1:out{i,j}.p, out{i,j}.cost,style1{k});
                    hold on;
                end
            end
        end
    end
end

if(0)
    for k=1:6
        eval(sprintf('out=out1%d;',k));
        rse = zeros(size(out));
        for i=1:size(out,1)
            for j=1:size(out,2)
                if(~isempty(out{i,j}))
                    figure((i-1)*size(out,2)+j);
                    semilogy(1:out{i,j}.p, out{i,j}.RMSE,style1{k});
                    hold on;
                end
            end
        end
    end
end

if(1)
    %% show the comparison between different method on different number of projections
    prj = zeros(6,1);
    intval = 6:-1:1;
    for i=1:length(intval)
        prj(i) = length((0:intval(i):179));
    end
    figure;
    for k=1:6
        eval(sprintf('out=out1%d;',k));
        rse = zeros(size(out));
        for i=1:size(out,1)
            for j=1:size(out,2)
                if(~isempty(out{i,j}))
                    rse(i,j)=min(out{i,j}.RMSE);
                end
            end
        end
        if(length(prj)==length(rse))
            semilogy(prj,rse,style2{k}); hold on;
        end
        eval(sprintf('rse2%d=rse;',k));
    end
    out=out1;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:<');

    out=out3;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'g:^');

    out=out5;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'c-.p');

    out=out2;
    rse = zeros(size(out));
    for i=1:size(out,1)
        for j=1:size(out,2)
            if(~isempty(out{i,j}))
                rse(i,j)=min(out{i,j}.RMSE);
            end
        end
    end
    semilogy(prj,min(rse'),'k--d'); hold on;
    eval(sprintf('rsel%d=rse;',k));

    out=out4;
    rse = ones(size(out));
    for i=1:size(out,1)
        for j=1:size(out,2)
            if(~isempty(out{i,j}))
                rse(i,j)=min(out{i,j}.RMSE);
            end
        end
    end
    semilogy(prj,min(rse'),'r--d'); hold on;
    eval(sprintf('rse%d=rse;',4));
end


if(0)
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

    %% Find out that when a=-6.5, it give the lowest RSE
    for k=1:6
        eval(sprintf('out=out2%d;',k));
        rse = zeros(size(out));
        for i=1:size(out,1)
            for j=1:size(out,2)
                if(~isempty(out{i,j}))
                    rse(i,j)=min(out{i,j}.RMSE);
                end
            end
        end
        eval(sprintf('rse2%d=rse;',k));
    end
end


