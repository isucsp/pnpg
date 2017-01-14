function compare(field,plotFunc,varargin)
    colors={'g','b','r','k','m'};
    marks = {'.','*','o','s','^'};
    lines = {'-','--','-.'};

    if(~iscell(field))
        fieldname{1}=field;
    else
        fieldname=field;
    end

    vars=varargin;

    m=0;
    figure;
    for ii=1:length(vars)
        if(ii==2)
            hold on;
        end
        for j=1:length(fieldname)
            v{j}=getfield(vars{ii},fieldname{j});
        end
        if(length(fieldname)==1)
            plotFunc(v{1}-m,[colors{mod(ii-1,length(colors))+1} lines{mod(ii-1,length(lines))+1}]);
        elseif(length(fieldname)==2)
            plotFunc(v{1},v{2}-m,[colors{mod(ii-1,length(colors))+1} lines{mod(ii-1,length(lines))+1}]);
            xlabel(fieldname{1});
            ylabel(fieldname{2});
        end
        if(isfield(vars{ii},'name'))
            names{ii}=filterName(getfield(vars{ii},'name'));
        else
            names{ii}=sprintf('%c',ii+96);
        end
    end
    if(length(fieldname)==1)
        title(fieldname{1});
    elseif(length(fieldname)==2)
        title([fieldname{1} ' v.s. ' fieldname{2}]);
    end
    legend(names);
end

function o = filterName(name)
    o=strrep(name,'_','\_');
end


