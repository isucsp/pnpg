function [data,dist,img]=readRawData(filename,slides,thresh,opt)
% read the "slides"-th slide from the raw data *.sin and *.raw
% thresh is provided to specify the how stric the boundary of the sinogram
% will be used for estimation of the rotation center and X-ray source to
% rotation center distance.
%
% One example to call is 
%     
%     readRawData('../../cnde/blade_a/blade%03d.sin',220);
% or  readRawData('../../cnde/blade_a/blade%03d.raw',220);


    [a,b]=regexpi(filename,'\.raw$');
    if(~isempty(a)) % read raw data
        str='';
        for i=0:359
            fi=fopen(sprintf(filename,i),'r');
            temp = reshape(fread(fi,inf,'uint16'),3072,[]);
            if(~isempty(temp))
                data(:,i+1)=temp(:,slides);
                temp=length(str);
                str=sprintf('reading the %d-th projection',i);
                fprintf([repmat('\b',1,temp) '%s'],str);
                % figure(1); showImg(-log(temp)); title(sprintf('%03d',i)); drawnow;
            end;
            fclose(fi);
        end
    end

    [a,b]=regexpi(filename,'\.sin$');
    if(~isempty(a)) % read sin data
        HEADERLINES = 5;
        DELIMITER = '\t';
        filename=sprintf(filename,floor(slides/10));
        file=fopen(filename,'r');
        for i=1:HEADERLINES
            fgetl(file);
        end
        temp=mod(slides,10)*360;
        while(temp>0)
            if(~isempty(strtrim(fgetl(file)))) temp=temp-1; end
        end

        i=1;
        while(i<=360)
            temp=strtrim(fgetl(file));
            if(~isempty(temp)) data(:,i)=str2num(temp)'; i=i+1; end;
        end
    end

    data(data==0)=min(data(data>0));
    I_max=max(data(:));
    data=-log(data/I_max);

    if(exist('thresh','var') && ~isempty(thresh))
        [data,dist,img,center]=preprocessCTdata(data,thresh);
    else
        [data,dist,img,center]=preprocessCTdata(data);
    end

    %imwrite(I/max(I(:)),sprintf('%s%02d%01d.jpg',prefix,i,j),'jpg');
    %save(sprintf('%s%02d%01d.mat',prefix,i,j),'I');
end

