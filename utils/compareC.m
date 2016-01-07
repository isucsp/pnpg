% centeralized compare
function compareC(field,plotFunc,varargin)
    colors={'g','b','r','k','m'};
    marks = {'.','*','o','s','^'};
    lines = {'-','--','-.'};

    if(~iscell(field))
        fieldname{1}=field;
    else
        fieldname=field;
    end

    m=+inf;
    for ii=1:length(varargin)
        for j=1:length(fieldname)
            v{j}=getfield(varargin{ii},fieldname{j});
            if(j==length(fieldname))
                m=min(m,min(v{j}));
            end
        end
    end
    figure;
    for ii=1:length(varargin)
        if(ii==2)
            hold on;
        end
        for j=1:length(fieldname)
            v{j}=getfield(varargin{ii},fieldname{j});
        end
        if(length(fieldname)==1)
            plotFunc(v{1}-m,[colors{mod(ii-1,length(colors))+1} lines{mod(ii-1,length(lines))+1}]);
        elseif(length(fieldname)==2)
            plotFunc(v{1},v{2}-m,[colors{mod(ii-1,length(colors))+1} lines{mod(ii-1,length(lines))+1}]);
            xlabel(fieldname{1});
            ylabel(fieldname{2});
        end
        str{ii}=sprintf('%c',ii+96);
    end
    if(length(fieldname)==1)
        title(['Centralized ' fieldname{1}]);
    elseif(length(fieldname)==2)
        title([fieldname{1} ' v.s. Centeralized ' fieldname{2}]);
    end
    legend(str);
end

