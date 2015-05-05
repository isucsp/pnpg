function [data,dist,img]=readRawData(filename,slides,thresh)
% read the "slides"-th slide from the raw data *.sin 
% thresh is provided to specify the how stric the boundary of the sinogram
% will be used for estimation of the rotation center and X-ray source to
% rotation center distance.
%
% One example to call is 
%     
%     readRawData('../../cnde/blade_a/blade1a%03d.sin',220);
%
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

    data(data==0)=min(data(data>0));
    I_max=max(data(:));
    data=-log(data/I_max);

    if(exist('thresh','var'))
        [data,dist,img,center]=preprocessCTdata(data,thresh);
    else
        [data,dist,img,center]=preprocessCTdata(data);
    end

    %imwrite(I/max(I(:)),sprintf('%s%02d%01d.jpg',prefix,i,j),'jpg');
    %save(sprintf('%s%02d%01d.mat',prefix,i,j),'I');
end
