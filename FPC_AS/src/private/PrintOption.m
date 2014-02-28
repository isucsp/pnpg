function PrintOption(opts)
%-------------------------------------------------------------------------
% print fields of "opts"
%-------------------------------------------------------------------------

name = fieldnames(opts); nstruct = size(name, 1);  
fprintf('\nOptions are: \n');
for indi = 1 : nstruct
    switch class(opts.(name{indi}))
        case 'char'
            fprintf(1, '%20s : %s',name{indi}, opts.(name{indi}));
            if isfield(opts.memo, name{indi});
                fprintf(1, '\t  %s',opts.memo.(name{indi}));
            end
            fprintf('\n');
            
        case 'double'

            temp = opts.(name{indi});
            [nx, ny] = size(temp);
            if nx > 1 || ny > 1
                fprintf(1, '%20s, (first 5 elements): ',name{indi});
                nx = min(nx, 5);     ny = min(ny, 5);  
                fprintf(1, '%g \t', temp(1:nx,1:ny)); 
            elseif nx == 1 && ny == 1
                fprintf(1, '%20s : %g',name{indi}, opts.(name{indi}));
            else
                fprintf(1, '%20s : \t',name{indi});
            end
            
            if isfield(opts.memo, name{indi});
                fprintf(1, '\t  %s',opts.memo.(name{indi}));
            end
            fprintf('\n');
    end

end
end
