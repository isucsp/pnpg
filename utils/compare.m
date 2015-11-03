
function compare(field,plotFunc,varargin)
    colors={'g','b','r','k','m'};
    marks = {'.','*','o','s','^'};
    lines = {'-','--','-.'};

    if(~iscell(field))
        fieldname{1}=field;
    else
        fieldname=field;
    end

    figure;
    m=+inf;
    for ii=1:length(varargin)
        if(ii==2)
            hold on;
        end
        for j=1:length(fieldname)
            v{j}=getfield(varargin{ii},fieldname{j});
            m=min(m,min(v{j}));
        end
        m=m-1e-2;
        if(length(fieldname)==1)
            plotFunc(v{1}-m,[colors{mod(ii-1,length(colors))+1} lines{mod(ii-1,length(lines))+1}]);
        elseif(length(fieldname)==2)
            plotFunc(v{1},v{2},[colors{mod(ii-1,length(colors))+1} lines{mod(ii-1,length(lines))+1}]);
        end
        str{ii}=sprintf('%d',ii);
    end
    if(length(fieldname)==1)
        title(fieldname{1});
    elseif(length(fieldname)==2)
        title([fieldname{1} ' v.s. ' fieldname{2}]);
    end
    legend(str);
end

