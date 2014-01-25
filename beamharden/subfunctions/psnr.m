function y = psnr(im1,im2,scale);

if nargin<3
    scale=1;
end

[m,n] = size(im1);
x1 = double(im1(:));
x2 = double(im2(:));
mse = norm(x1-x2);
mse = (mse*mse)/(m*n);
if mse >0
    temp=double(scale^2/mse);
    y = 10*log10(temp);
else
    disp('infinite psnr');
    y=1e30;
end
