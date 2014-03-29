function summary = plotVector(out,oMode,field,plotf)
    color  = {'r' ,'g' ,'b' ,'k' ,'c' ,'y' ,'m'};
    lines  = {'-' ,'--',':' ,'-.'};
    style2 = {'b-*','r-*','b:o','r:o','b--s','r--s'};

    if(oMode==0)
        maxlen=1;
        figure;
        for i=1:length(out(:))
            if(~isempty(out{i}))
                res{i}=getfield(out{i},field);
                plotf(res{i},'color',color{mod(i-1,length(color))+1},...
                    'linestyle',lines{mod(i-1,length(lines))+1});
                hold on;
                temp=res{i};

                if(length(temp(:))>maxlen) maxlen=length(temp(:)); end
            else
                plotf(0); hold on;
            end
            str{i}=num2str(i);
        end
        legend(str);
        summary=zeros(maxlen, length(out(:)));
        for i=1:length(out(:))
            temp=res{i};
            summary(1:length(temp),i)=temp(:);
        end
    elseif(oMode==1)
        for i=1:length(out(:))
            if(~isempty(out{i}))
                res{i}=getfield(out{i},field);
                plotf(res{i},'color',color{mod(i-1,length(color))+1},...
                    'linestyle',lines{mod(i-1,length(lines))+1});
                hold on;
                temp=res{i};

                if(length(temp(:))>maxlen) maxlen=length(temp(:)); end
            else
                plotf(0); hold on;
            end
            str{i}=num2str(i);
        end

    end

    
end
