function img = remoteMatlab(n)
    !scp co3125-07:research/phi/src/reImg.data ./
    f=fopen('reImg.data','r');
    img=reshape(fread(f,n^2,'float'),n,n)';
    figure; imshow(img,[]);
end

