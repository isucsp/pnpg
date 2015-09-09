function img = localMatlab(fn,col)
        f=fopen(fn,'r');
        img=fread(f,inf,'float');
        if(exist('col','var'))
            img=reshape(img,[],col);
        end
        figure; showImg(img);
end

