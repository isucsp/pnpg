function [newdata,dist,img,center]=preprocessCTdata(data,thresh)
% Syntax:
%     [newdata,dist,img,center]=preprocessCTdata(data,thresh)
%
% this subroutine is used to find the center of the projection and the
% distance from the rotation center to X-ray source.
% Here, "data" is after taking the logarithm
%
% thresh is provided to specify the how stric the boundary of the sinogram
% will be used for estimation of the rotation center and X-ray source to
% rotation center distance.
%
% thresh is used in for the sinogram after taking the logarithm with value
% in (0,1), where 0 and 1 correspond to the minimum and maximum of the
% sinogram.
%
% Outputs:
% dist:        estimated distance from the X-ray source to the rotation center
% img:         one fan-beam FBP reconstruction
% center:      the estimated position of the projection of rotation center on the detector array.
%              (this is not very useful here, since newdata is already
%              adjusted and centerized accordingly.)
%
% author: Renliang Gu (gurenliang@gmail.com)
%

if(nargin==0)
    help preprocessCTdata
end

if(~exist('thresh','var'))
    thresh=0.2;
    fprintf('No input for thresh, use %g as default\n',thresh);
end;

figs=[];

while(true)

    thresh=thresh*max(data(:));

    mask=data>thresh;
    [N,M]=size(mask);
    for i=1:M
        g1(i,1) = find(mask(:,i),1);
        g2(i,1) = find(mask(:,i),1,'last');
    end
    center=(mean(g1)+mean(g2))/2;
    figs=[figs; figure];
    showImg(data); hold on; plot(g1,'r'); plot(g2,'r');
    drawnow;

    gl=g2-center; g2=g1-center; g1=gl; clear 'gl';
    G1=fft(g1); G2=fft(g2); gg=ifft(conj(G1).*(G2));
    dd(1)=mean(g1)/tan((find(gg==min(gg),1)-M/2)/M*2*pi);
    dd(2)=max(g1)/tan(...
        mod(find(g2==min(g2),1)-find(g1==max(g1),1)+M/2,M)*2*pi/M );
    dd(3)=min(g1)/tan(...
        mod(find(g2==max(g2),1)-find(g1==min(g1),1)+M/2,M)*2*pi/M );

    theta=(0:M-1)*2*pi/M; theta=theta(:);
    gg2=@(dddd) interp1([theta; theta+2*pi], [g2; g2],...
        theta+2*atan(  g1/dddd  )+pi,'spline');
    objfunc=@(ddd) norm(g1+gg2(ddd));
    [dist,~,status]=fminsearch(objfunc,median(dd));
    if(status~=1) keyboard; end

    figs=[figs; figure];
    fplot(objfunc,dist+[-100,100]); hold on; plot(dist,objfunc(dist),'r*');

    figs=[figs; figure];
    subplot(2,1,1); plot(theta,g1,'r'); hold on;
    plot(theta,-gg2(dist),'b');
    subplot(2,1,2); plot(theta,g1+gg2(dist));

    band=floor(min(min(center-2,N-1-center),max(g1)*1.2));
    newdata=[];
    for i=1:M
        newdata(:,i) = interp1((1:N), data(:,i), (center-band):(center+band+1),'linear');
    end

    img=ifanbeam(newdata(end:-1:1,:),dist,...
        'FanCoverage','cycle','FanRotationIncrement',1,...
        'FanSensorGeometry','line','FanSensorSpacing',1,...
        'OutputSize',length(newdata(:,1)));
    figs=[figs; figure];
    showImg(img,0);

    while(true)
        fprintf('Estimated distance from X-ray source to rotation center is %g times pixel size\n', dist);
        str=sprintf('enter 0 to accept the current thresh=%g,\nenter 1 to change value of thresh\nenter 2 to enter debug mode:',thresh/max(data(:)));
        decide=input(str);
        switch(decide)
            case 0
                return;
            case 1
                thresh=input('enter new desired value for thresh:');
                break;
            case 2
                fprintf('enter ''dbcont'' or press ''F5'' to continue\n');
                keyboard
        end
    end
end

end

