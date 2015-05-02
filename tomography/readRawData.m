function [data,dist,img]=readRawData(filename,slides,thresh)
    DELIMITER = '\t';
    HEADERLINES = 5;
    % Import the file
    newData=importdata(sprintf(filename,floor(slides/10)),DELIMITER,HEADERLINES);

    j=mod(slides,10);
    data=newData.data';
    data=data(:,j*360+1:(j+1)*360);
    data(data==0)=min(data(data>0));
    I_max=max(data(:));
    data=-log(data/I_max);

    [data,dist,img,center]=preprocessCTdata(data);

    %imwrite(I/max(I(:)),sprintf('%s%02d%01d.jpg',prefix,i,j),'jpg');
    %save(sprintf('%s%02d%01d.mat',prefix,i,j),'I');
end
