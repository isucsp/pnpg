function splitmat(filename,opt)
    % rescue the corrupted mat file and save its variables one be one into 
    % splitted mat files.
    % if opt=0 or unspecified, then find the variables until the corrupted
    % one
    % if opt=1, find possible variable with type 0x0F 00 00 00
     
    if(~exist('opt','var')) opt=0; end

    [a,~]=regexpi(filename,'\.mat$');
    if(~isempty(a)) % read raw data
        filename=filename(1:a-1);
    end
    f=fopen([filename '.mat'], 'rb');
    header=fread(f,128);
    count=0;

    switch(opt)
        case 0
            while true
                count = count+1;
                h2=fread(f,2,'int32');
                if length(h2) < 2
                    disp 'finished reading file';
                    break;
                end
                if h2(2) == 0
                    fprintf('h2=%s\n',num2str(h2));
                    fprintf('Found bad 0-byte size at variable #%d.', count);
                    break;
                end
                fout=fopen(sprintf('%s_%02d.mat', filename, count), 'wb');
                fwrite(fout, header);
                fwrite(fout, h2, 'int32');
                data = fread(f, h2(2));
                fwrite(fout, data);
                fclose(fout);
            end
        case 1
            h2=fread(f);
            idx=find(h2==15);
            i=0;
            strlen=0;

            while(i<length(idx))
                i=i+1;
                str=sprintf('%02d%%... ',floor(100*i/length(idx)));
                if((idx(i)+7)<=length(h2) && all(h2(idx(i)+1:idx(i)+3)==0))
                    fseek(f,128+idx(i)-1,-1);

                    h3=fread(f,2,'int32');
                    if length(h3) < 2
                        fprintf('i=%d, idx(i)=%d, h3=%s\n',i,idx(i),num2str(h3));
                        disp 'finished reading file';
                        strlen=0;
                        continue;
                    end

                    if h3(2) <= 0
                        % fprintf('i=%d/%d, idx(i)=%d, h3=%s\n',i,length(idx),idx(i),num2str(h3));
                        % disp(sprintf('Found bad 0-byte size at variable #%d.', i));
                        % strlen=0;
                        continue;
                    end

                    count=count+1;
                    newfilename=sprintf('%s_%02d.mat', filename, count);
                    fout=fopen(['/tmp/' newfilename], 'wb');
                    fwrite(fout, header);
                    fwrite(fout, h3, 'int32');
                    data = fread(f, h3(2));
                    fwrite(fout, data);
                    fclose(fout);
                    try
                        w=warning('off','all');
                        ss=load(['/tmp/' newfilename]);
                        warning(w);
                    catch
                        warning(w);
                        count=count-1;
                        continue;
                    end
                    i=find(idx<idx(i)+8+h3(2),1,'last');
                    system(['mv /tmp/' newfilename '  ./']);
                    str=sprintf('%s Recovered %s.\n',str,newfilename);
                end
                fprintf([repmat('\b',1,strlen) '%s'],str);
                if(str(end-2)=='t')
                    strlen=0;
                else
                    strlen = length(str);
                end
            end
    end
    fclose(f);
    fprintf('\n');

    temp___{1}=filename;
    temp___{2}=count;
    decide=input(sprintf('Combine all file to %s_rec.mat [y/N]?',filename),'s');
    if strcmpi(decide,'y')
        clear('-regexp','(?!temp___)^.*$');
        for i=1:temp___{2}
            load(sprintf('%s_%02d.mat', temp___{1}, i));
        end
        save(sprintf('%s_rec.mat', temp___{1}));
        for i=1:temp___{2}
            system(sprintf('rm %s_%02d.mat', temp___{1}, i));
        end
    end
end

