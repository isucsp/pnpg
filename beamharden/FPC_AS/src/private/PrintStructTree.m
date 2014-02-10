function PrintStructTree(structure,filename);


name = fieldnames(structure);

nstruct = size(name, 1);

fprintf(filename, '\n');

for indi = 1 : nstruct
    
    switch class(getfield(structure, name{indi}))
        case 'char'
            
%             name{indi}
%             class(getfield(structure, name{indi}))
%             getfield(structure, name{indi})
            
             fprintf(filename, '%30s : %s \n',name{indi}, getfield(structure, name{indi}));
%              fprintf(filename, '%15s :  \n',name{indi}); 
%              fprintf(filename, '%15s : %s \n',name{indi}, getfield(structure, name{indi}));

        case 'double'
%             name{indi}
%             class(getfield(structure, name{indi}))
%             getfield(structure, name{indi})
%             fprintf(filename, '%15s : %d \n',name{indi}, getfield(structure, name{indi}));
              
              temp = getfield(structure, name{indi});
              [nx, ny] = size(temp);
              if nx > 1 && ny > 1
                  fprintf(filename, '%30s :  \n',name{indi});
                  nx = min(nx, 50);     ny = min(ny, 50);
                  fprintf(filename, '%30d \n', temp(1:nx,1:ny));
              elseif nx == 1 && ny == 1
                  fprintf(filename, '%30s : %d \n',name{indi}, getfield(structure, name{indi}));    
              else
                  fprintf(filename, '%30s : \n',name{indi});
              end
              

    end
    
end

fprintf(filename,'\n');