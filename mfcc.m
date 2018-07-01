function [ CC, FBE, frames ] = mfcc( speech, fs, Tw, Ts, alpha, window, R, M, N, L )
% MFCC Mel frequency cepstral coefficient feature extraction.
%
%   MFCC(S,FS,TW,TS,ALPHA,WINDOW,R,M,N,L) returns mel frequency 
%   cepstral coefficients (MFCCs) computed from signal given 
%   in vector S and sampled at FS (Hz). The signal is first 
%   preemphasised using a first order FIR filter with preemphasis 
%   coefficient ALPHA. The preemphasised signal is subjected 
%   to the short-time Fourier transform analysis with frame durations 
%   of TW (ms), frame shifts of TS (ms) and analysis window function 
%   given as a function handle in WINDOW. This is followed by magnitude 
%   spectrum computation followed by filterbank design with M triangular 
%   filters uniformly spaced on the mel scale between lower and upper 
%   frequency limits given in R (Hz). The filterbank is applied to 
%   the magnitude spectrum values to produce filterbank energies (FBEs) 
%   (M per frame). Log-compressed FBEs are then decorrelated using the 
%   discrete cosine transform to produce cepstral coefficients. Final
%   step applies sinusoidal lifter to produce liftered MFCCs.
%
%   [CC,FBE,FRAMES]=MFCC(...) also returns FBEs and windowed frames,
%   with feature vectors and frames as columns.
%
    %% PRELIMINARIES 

    % Ensure correct number of inputs
    if( nargin~= 10 ), help mfcc; return; end; 

    % Explode samples to the range of 16 bit shorts
    if( max(abs(speech))<=1 ), speech = speech * 2^15; end;

    Nw = round( 1E-3*Tw*fs );    % frame duration (samples)
    Ns = round( 1E-3*Ts*fs );    % frame shift (samples)

    nfft = 2^nextpow2( Nw );     % length of FFT analysis 
    K = nfft/2+1;                % length of the unique part of the FFT 


    %% HANDY INLINE FUNCTION HANDLES

    % Forward and backward mel frequency warping 
    hz2mel = @( hz )( 1127*log(1+hz/700) );     % Hertz to mel warping function
    mel2hz = @( mel )( 700*exp(mel/1127)-700 ); % mel to Hertz warping function

    % Type III DCT matrix routine 
    dctm = @( N, M )( sqrt(2.0/M) * cos( repmat([0:N-1].',1,M) ...
                                       .* repmat(pi*([1:M]-0.5)/M,N,1) ) );

    % Cepstral lifter routine 
    ceplifter = @( N, L )( 1+0.5*L*sin(pi*[0:N-1]/L) );

    
    % Preemphasis filtering 
    speech = filter( [1 -alpha], 1, speech ); % fvtool( [1 -alpha], 1 );

    % Framing and windowing (frames as columns)
    frames = vec2frames( speech, Nw, Ns, 'cols', window, false );

    % Magnitude spectrum computation (as column vectors)
    MAG = abs( fft(frames,nfft,1) ); 

    % Triangular filterbank with uniformly spaced filters on mel scale
    H = trifbank( M, K, R, fs, hz2mel, mel2hz ); % size of H is M x K 

    % Filterbank application to unique part of the magnitude spectrum
    FBE = H * MAG(1:K,:); % FBE( FBE<1.0 ) = 1.0; % apply mel floor

    % DCT matrix computation
    DCT = dctm( N, M );

    % Conversion of logFBEs to cepstral coefficients through DCT
    CC =  DCT * log( FBE );

    % Cepstral lifter computation
    lifter = ceplifter( N, L );

    % Cepstral liftering gives liftered cepstral coefficients
    CC = diag( lifter ) * CC; % ~ HTK's MFCCs
end