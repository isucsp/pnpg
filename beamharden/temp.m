if(1)
r=size(out010,1);
c=size(out010,2);
for i=1:r
    for j=1:c
        if((i<=size(out,1) && j<=size(out,2) && ~isempty(out{i,j})))
            if((~isempty(out010{i,j}))) disp([i,j]); end
            continue;
        end
        out{i,j}=out010{i,j};
    end
end
end
