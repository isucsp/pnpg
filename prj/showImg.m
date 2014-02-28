function []=showImg(x, varargin)

while(~isempty(varargin) && iscell(varargin{1}) && length(varargin)==1)
    varargin=varargin{1};
end

if(length(varargin)==2)
    minI=varargin{1}; maxI=varargin{2};
elseif(length(varargin)==1)
    minI=varargin{1}; maxI=double(max(x(:)));
else
    minI=double(min(x(:))); maxI=double(max(x(:)));
end

if(isempty(minI)) minI=double(min(x(:))); end
if(isempty(maxI)) maxI=double(max(x(:))); end

if(min(size(x))==1) x=reshape(x,sqrt(length(x)),[]); end

imshow(x,[minI, maxI],'InitialMagnification','fit');
%imagesc(x);
colormap(gray)
axis off;
%set(gca,'position',[0,0,1,1]);
colorbar;
drawnow

end

