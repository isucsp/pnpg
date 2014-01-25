function img = localMatlab(fn)
        f=fopen(fn,'r');
        img=[];
        while(~feof(f))
            img=[img; fread(f,1000,'float')];
        end
        if(mod(length(img(:)),360)==0)
            img=reshape(img,[],360)';
        else
            img=reshape(img,sqrt(length(img(:))),[])';
        end
        figure; imshow(img,[]);
end

