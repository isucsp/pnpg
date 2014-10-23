%%%
%%%  Sub-Function:  designFilter
%%%

function filt = designFilter(filter, len, Ts)
% Returns the Fourier Transform of the filter which will be 
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter 
%               len    - the length of the projections
%
% OUTPUT ARGS:  filt   - the filter to use on the projections

order = len;
d=1; 

% First create a ramp filter - go up to the next highest
% power of 2.

filt = (0:2/order:1)';
w = pi*filt;   % frequency axis up to Nyquist 

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end).*(sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end).*cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end).*(.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end).*(1+cos(w(2:end)./d)) / 2;
    case 'renliang'
        N=order;
        n=-N/2:N/2-1;
        hn=-1./((pi*Ts*n').^2);
        hn(find(mod(n,2)==0))=0;
        hn(N/2+1)=1/4/Ts^2;
        hamming=0.54+0.46*cos(2*pi*n'/N);
        hann=0.5+0.5*cos(2*pi*n'/N);
        sinc=sin(pi*n'/N)./(pi*n'/N);
        sinc(N/2+1)=1;
        hn=fftshift(fft(fftshift(hn)))*Ts;
        hn=hn.*hann;
        hn=fftshift(ifft(fftshift(hn)))/Ts;
        hn=real(hn);
        hn=[zeros(N/2,1); hn; zeros(N/2,1)];
        N=N*2;
        Hf=fft(fftshift(hn))*Ts;
        filt=real(Hf(1:N/2+1));
    case 'renliang1'
        N=order;
        n=-N/2:N/2-1;
        hn=-1./((pi*Ts*n').^2);
        hn(find(mod(n,2)==0))=0;
        hn(N/2+1)=1/4/Ts^2;
        hamming=0.54+0.46*cos(2*pi*n'/N);
        hann=0.5+0.5*cos(2*pi*n'/N);
        sinc=sin(pi*n'/N)./(pi*n'/N);
        sinc(N/2+1)=1;
        hn=fftshift(fft(fftshift(hn)))*Ts;
        hn=hn.*hann;
        filt=real(hn(1:N/2+1));
    otherwise
        eid = sprintf('Images:%s:invalidFilter',mfilename);
        msg = 'Invalid filter selected.';
        error(eid,'%s',msg);
end

filt = [filt(:) ; filt(end-1:-1:2)];    % Symmetry of the filter
