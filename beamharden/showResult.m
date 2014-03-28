function res = plotVector(out,field,plotf)
    color  = {'r' ,'g' ,'b' ,'k' ,'c' ,'y' ,'m'};
    lines  = {'-' ,'--',':' ,'-.'};
    style2 = {'b-*','r-*','b:o','r:o','b--s','r--s'};

    figure;
    for i=1:length(out)
        res{i}=getfield(out{i},field);
        plotf(res{i},'color',color{mod(i,length(color))+1},...
            'linestyle',lines{mod(i,length(lines))+1});
        hold on;
        str{i}=num2str(i);
    end
    legend(str);
end
