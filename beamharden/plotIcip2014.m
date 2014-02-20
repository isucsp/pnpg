
style1 = {'b-','r-','b:','r:','b--','r--','-','*','.','+','d'};
style2 = {'b-*','r-*','b:o','r:o','b--s','r--s'};

%% NDE report plot
if(1)
    aa=0;
    prj = zeros(6,1);
    intval = 6:-1:1;
    for i=1:length(intval)
        prj(i) = length((0:intval(i):179));
    end
    figure;
    for k=1:2:6
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
            aa=aa+1;
            legStr{aa} = ['unknown Ie by ' num2str(k)];
        end
        eval(sprintf('rse2%d=rse;',k));
    end
    out=out1;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:<');
    aa=aa+1; legStr{aa}='known Ie by dis';
    out=out6;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:>','linewidth',2);
    aa=aa+1; legStr{aa}='known Ie by b0';
    out=out7;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:h');
    aa=aa+1; legStr{aa}='known Ie by b1';

    out=out3;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'g:^');
    aa=aa+1; legStr{aa}='direct FBP';

    out=out5;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'c-.p');
    aa=aa+1; legStr{aa}='FBP after linearization';

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
    aa=aa+1; legStr{aa}='FPCAS';

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
    aa=aa+1; legStr{aa}='FPCAS after linearization';

    legend(legStr);
    for k=1:2:6
        eval(sprintf('out=out1%d;',k));
        rse = zeros(size(out));
        for i=1:size(out,1)
            for j=1:size(out,2)
                if(~isempty(out{i,j}))
                    figure((i-1)*size(out,2)+j);
                    subplot(2,1,1);
                    loglog(1:out{i,j}.p, out{i,j}.cost,style1{k}); hold on;
                    subplot(2,1,2);
                    loglog(1:out{i,j}.p, out{i,j}.RMSE,style1{k});
                    hold on;
                end
            end
        end
    end
end

if(0)
    %% show the comparison between different method on different number of projections
    aa=0;
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
            aa=aa+1;
            legStr{aa} = ['unknown Ie by ' num2str(k)];
        end
        eval(sprintf('rse2%d=rse;',k));
    end
    out=out1;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:<');
    aa=aa+1; legStr{aa}='known Ie by dis';
    out=out6;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:>','linewidth',2);
    aa=aa+1; legStr{aa}='known Ie by b0';
    out=out7;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'r:h');
    aa=aa+1; legStr{aa}='known Ie by b1';

    out=out3;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'g:^');
    aa=aa+1; legStr{aa}='direct FBP';

    out=out5;
    rse=zeros(size(out));
    for i=1:length(out)
        rse(i)=min(out{i}.RMSE);
    end
    plot(prj,rse,'c-.p');
    aa=aa+1; legStr{aa}='FBP after linearization';

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
    aa=aa+1; legStr{aa}='FPCAS';

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
    aa=aa+1; legStr{aa}='FPCAS after linearization';

    legend(legStr);
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

%% runlist = 30
if(0)
    figure;
    out=out30; p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.RMSE,style1{k});
        p=p+out{k}.p;
        hold on;
    end
end

%% runlist = 31
if(0)
    figure;
    out=out31; p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.RMSE,style1{1});
        p=p+out{k}.p;
        hold on;
    end
    p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.l1Pen,style1{2});
        p=p+out{k}.p;
        hold on;
    end
    p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.cost,style1{3});
        p=p+out{k}.p;
        hold on;
    end
    p=1;
    for k=1:length(out)
        semilogy(p:p+out{k}.p-1, out{k}.nonneg,style1{4});
        p=p+out{k}.p;
        hold on;
    end
end

