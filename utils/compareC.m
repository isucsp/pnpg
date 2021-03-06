% centeralized compare
function compareC(field,plotFunc,varargin)
    colors= {'g','b','r','k','m'};
    marks = {'+','x','d','v','<','>','p','h','.','*','o','s','^'};
    lines = {'-','--','-.',':'};

    if(~iscell(field))
        fieldname{1}=field;
    else
        fieldname=field;
    end

    vars=varargin;

    m=+inf;
    for ii=1:length(vars)
        m=min(m,min(vars{ii}.(fieldname{end})));
    end
    figure;
    for ii=1:length(vars)
        if(ii==2)
            hold on;
        end
        for j=1:length(fieldname)
            v{j}=vars{ii}.(fieldname{j});
        end
        colorIdx=mod(ii-1,length(colors))+1;
        tmp=floor((ii-1)/length(colors))+1;
        lineIdx=mod(tmp-1,length(lines))+1;
        if(length(fieldname)==1)
            plotFunc((v{1}-m)/m,[colors{colorIdx} lines{lineIdx}]);
            ylabel(fieldname{1});
        elseif(length(fieldname)==2)
            plotFunc(v{1},(v{2}-m)/m,[colors{colorIdx} lines{lineIdx}]);
            xlabel(fieldname{1});
            ylabel(fieldname{2});
        end
        if(isfield(vars{ii},'name'))
            names{ii}=filterName(vars{ii}.('name'));
        else
            names{ii}=sprintf('%c',ii+96);
        end
    end
    if(length(fieldname)==1)
        title(['Centralized ' fieldname{1}]);
    elseif(length(fieldname)==2)
        title([fieldname{1} ' v.s. Centeralized ' fieldname{2}]);
    end
    legend(names);
end

function o = filterName(name)
    o=strrep(name,'_','\_');
end


