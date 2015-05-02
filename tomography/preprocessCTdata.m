
function [newdata,dist,img,center]=preprocessCTdata(data,theta)
    % this subroutine is used to find the center of the projection and the
    % distance from the rotation center to X-ray source.
    
    if(~exist('thresh','var')) thresh=0.5; end;

    mask=data>thresh;
    [N,M]=size(mask);
    for i=1:M
        g1(i,1) = find(mask(:,i),1);
        g2(i,1) = find(mask(:,i),1,'last');
    end
    center=(mean(g1)+mean(g2))/2;
    figure; showImg(data);
    figure; showImg(data); hold on; plot(g1,'r'); plot(g2,'r');

    gl=g2-center; g2=g1-center; g1=gl; clear 'gl';
    G1=fft(g1); G2=fft(g2); gg=ifft(G1.*conj(G2));
    dd(1)=mean(g1)/tan((find(gg==min(gg))-M/2)/M*2*pi);
    dd(2)=max(g1)/tan((find(g1==max(g1),1)-find(g2==min(g2),1)-M/2)*2*pi/M);
    dd(3)=min(g1)/tan((find(g1==min(g1),1)-find(g2==max(g2),1)-M/2)*2*pi/M);
    ddRange=floor(min(abs(dd))):ceil(max(abs(dd)));

    theta=(0:M-1)*2*pi/M; theta=theta(:);
    gg2=@(dddd) interp1([theta; theta+2*pi], [g2; g2],...
        theta+2*atan(  g1/dddd  )+pi,'spline');
    for i=1:length(ddRange)
        cost(i)=norm(g1+gg2(ddRange(i)));
    end
    dist=ddRange(find(cost==min(cost),1));

    figure; subplot(2,1,1); plot(theta,g1,'r'); hold on;
    plot(theta,-gg2(dist),'b');
    subplot(2,1,2); plot(theta,g1+gg2(dist));

    band=min(center-2,N-1-center);
    for i=1:M
        newdata(:,i) = interp1((1:N), data(:,i), (center-band):(center+band+1),'linear');
    end

    img=ifanbeam(newdata,dist,...
        'FanCoverage','cycle','FanRotationIncrement',1,...
        'FanSensorGeometry','line','FanSensorSpacing',1,...
        'OutputSize',length(newdata(:,1)));
    figure; showImg(img,0);

end

