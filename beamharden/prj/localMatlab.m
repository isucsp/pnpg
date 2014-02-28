function img = localMatlab(fn,col)
        f=fopen(fn,'r');
        img=[];
        while(~feof(f))
            img=[img; fread(f,1000,'float')];
        end
        img=reshape(img,[],col);
        figure; imshow(img,[]);
end

