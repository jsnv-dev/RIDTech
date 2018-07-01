function [ H, f, c ] = trifbank( M, K, R, fs, h2w, w2h )
% TRIFBANK Triangular filterbank.
%
%   [H,F,C]=TRIFBANK(M,K,R,FS,H2W,W2H) returns matrix of M triangular filters 
%   (one per row), each K coefficients long along with a K coefficient long 
%   frequency vector F and M+2 coefficient long cutoff frequency vector C. 
%   The triangular filters are between limits given in R (Hz) and are 
%   uniformly spaced on a warped scale defined by forward (H2W) and backward 
%   (W2H) warping functions.
%  
    if( nargin~= 6 ), help trifbank; return; end; % very lite input validation

    f_min = 0;          % filter coefficients start at this frequency (Hz)
    f_low = R(1);       % lower cutoff frequency (Hz) for the filterbank 
    f_high = R(2);      % upper cutoff frequency (Hz) for the filterbank 
    f_max = 0.5*fs;     % filter coefficients end at this frequency (Hz)
    f = linspace( f_min, f_max, K ); % frequency range (Hz), size 1xK
    fw = h2w( f );

    % filter cutoff frequencies (Hz) for all filters, size 1x(M+2)
    c = w2h( h2w(f_low)+[0:M+1]*((h2w(f_high)-h2w(f_low))/(M+1)) );
    cw = h2w( c );

    H = zeros( M, K );                  % zero otherwise
    for m = 1:M 

        % implements Eq. (6.140) on page 314 of [1] 
        % k = f>=c(m)&f<=c(m+1); % up-slope
        % H(m,k) = 2*(f(k)-c(m)) / ((c(m+2)-c(m))*(c(m+1)-c(m)));
        % k = f>=c(m+1)&f<=c(m+2); % down-slope
        % H(m,k) = 2*(c(m+2)-f(k)) / ((c(m+2)-c(m))*(c(m+2)-c(m+1)));

        % implements Eq. (6.141) on page 315 of [1]
        k = f>=c(m)&f<=c(m+1); % up-slope
        H(m,k) = (f(k)-c(m))/(c(m+1)-c(m));
        k = f>=c(m+1)&f<=c(m+2); % down-slope
        H(m,k) = (c(m+2)-f(k))/(c(m+2)-c(m+1));
       
   end

   % H = H./repmat(max(H,[],2),1,K);  % normalize to unit height (inherently done)
   % H = H./repmat(trapz(f,H,2),1,K); % normalize to unit area 
end