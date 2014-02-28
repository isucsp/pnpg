function out = showImgMask(img,mask,varargin)
    Img=zeros(size(mask));
    maskIdx=find(mask);
    if(length(maskIdx)~=length(img(:)))
        Img=reshape(img,sqrt(length(img(:))),[]);
    else
        Img(maskIdx)=img;
    end
    if(nargout==0)
        showImg(Img,varargin);
    else
        out=Img;
    end
end
