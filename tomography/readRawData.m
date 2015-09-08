function [sino,dist,img]=readRawData(filename,slide)
% Syntax:
%
% [sino,dist,img]=readRawData(filename,slide)
%
% read the "slide"-th slide from the raw data *.sin and *.raw
%
% One example to call is 
%     
%     readRawData('../../cnde/blade_a/blade%03d.sin',220);
% or  readRawData('../../cnde/blade_a/blade%03d.raw',220);
%
% It returns the sinogram data after the logarithm (sino), the distance from the
% X-ray source to the rotation center (dist), and one FBP reconstrction of
% the "slide"-th slide of the data.
%
% Author: Renliang Gu (gurenliang@gmail.com)
%

if(nargin==0)
    help readRawData
end

file=[];

[a,b]=regexpi(filename,'\.raw$');
if(~isempty(a)) % read raw data
    prjWidth=3072;  % here we suppose the width of detector is 3072
    str='';
    for i=0:359
        file=fopen(sprintf(filename,i),'r');
        if(0==fseek(file,prjWidth*(slide-1)*2,'bof'))
            temp=fread(file,prjWidth,'uint16');
            if(~isempty(temp))
                Imea(:,i+1)=temp(:);
                temp=length(str);
                str=sprintf('reading the %d-th projection',i);
                fprintf([repmat('\b',1,temp) '%s'],str);
                % figure(1); showImg(-log(temp)); title(sprintf('%03d',i)); drawnow;
            else
                fprintf('error while reading %s\n',sprintf(filename,i));
                return;
            end;
        end
        fclose(file);
    end
    fprintf('\n');
end

[a,b]=regexpi(filename,'\.sin$');
if(~isempty(a)) % read sin data
    HEADERLINES = 5;
    DELIMITER = '\t';
    file=fopen(sprintf(filename,floor(slide/10)),'r');
    for i=1:HEADERLINES
        fgetl(file);
    end
    temp=mod(slide,10)*360;
    while(temp>0)
        if(~isempty(strtrim(fgetl(file)))) temp=temp-1; end
    end

    i=1;
    while(i<=360)
        temp=strtrim(fgetl(file));
        if(~isempty(temp)) Imea(:,i)=str2num(temp)'; i=i+1; end;
    end
    fclose(file);
end

if(isempty(file))
    fprintf('invalid filename!\n Return!\n');
    return;
end

Imea(Imea==0)=min(Imea(Imea>0));
I_max=max(Imea(:));
sino=-log(Imea/I_max);

[sino,dist,img,center]=preprocessCTdata(sino);

end

