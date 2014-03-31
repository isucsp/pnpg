if(1)
r=size(out012,1);
c=size(out012,2);
for i=1:r
    for j=1:c
        if((i<=size(out,1) && j<=size(out,2) && ~isempty(out{i,j})))
            if((~isempty(out012{i,j}))) disp([i,j]); end
            continue;
        end
        out{i,j}=out012{i,j};
    end
end
end
