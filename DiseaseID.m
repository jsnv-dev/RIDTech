function microbes = DiseaseID(AmbienceFileFullPath, MicrobesFileFullPath);
%% DISEASE IDENTIFICATION
%
% This program identifies microoragnisms that may cause rice plant diseases 
% (Magnaporthe oryzae, Thanatephorus cucumeris,and Xanthomonas oryzae which
% cause Leaf Blast, Sheath Blight, and Bacterial Leaf Blight respectively)based on 
% mel freqeuncy cepstral coefficients(MFCC) using sound signal processing 
% techniques and pattern recognition fuzzy neural network.

Tw = 25;                              % analysis frame duration (ms)
Ts = 10;                              % analysis frame shift (ms)
alpha = 0.97;                         % preemphasis coefficient
M = 20;                               % number of filterbank channels 
C = 12;                               % number of cepstral coefficients
L = 22;                               % cepstral sine lifter parameter
LF = 300;                             % lower frequency limit (Hz)
HF = 100000;                          % upper frequency limit (Hz)
fs = 200000;                          % sampling frequency
n = 24;                               % no. of bits
x1 = audioread(AmbienceFileFullPath); % file path of the ambient noise
x2 = audioread(MicrobesFileFullPath); % file path of the sound with microbes
x1 = x1(650001:end);
x2 = x2(400001:600000);
signal = [x1;x2];

%% Spectral Subtraction
%
% Short time Spectral Amplitude Minimum Mean Square Error Method for 
% Denoising noisy signal. 'signal' is the input noisy signal,
% IS is the initial silence (The initial time of the signal is 
% inactive period and may be used for initial noise parameter.)
% The output is the restored clean signal. 
%
 
if (nargin<3 | isstruct(IS)) 
    IS=.25;     % Initial Silence or Noise Only part in seconds 
end 
W=fix(.025*fs); % Window length is 25 ms 
SP=.4;          % Shift percentage is 40% (10ms) %Overlap-Add method works good with this value(.4) 
wnd=hamming(W); 
 
if (nargin>=3 & isstruct(IS))
    W=IS.windowsize 
    SP=IS.shiftsize/W; 
    %nfft=IS.nfft; 
    wnd=IS.window; 
    if isfield(IS,'IS') 
        IS=IS.IS; 
    else 
        IS=.25; 
    end 
end 

pre_emph=0; 
signal=filter([1 -pre_emph],1,signal); 
 
NIS=fix((IS*fs-W)/(SP*W) +1);       % number of initial silence segments 
y=segment(signal,W,SP,wnd);         % This function chops the signal into frames 
Y=fft(y); 
YPhase=angle(Y(1:fix(end/2)+1,:));  % Noisy Signal Phase 
Y=abs(Y(1:fix(end/2)+1,:));         % Spectogram 
numberOfFrames=size(Y,2); 
FreqResol=size(Y,1); 
 
N=mean(Y(:,1:NIS)')';               % initial Noise Power Spectrum mean 
LambdaD=mean((Y(:,1:NIS)').^2)';    % initial Noise Power Spectrum variance 
alpha=.99;                          % used in smoothing xi (For Deciesion Directed method for estimation of A Priori SNR) 
 
NoiseCounter=0; 
NoiseLength=9;                      % This is a smoothing factor for the noise updating 
 
G=ones(size(N));                    % Initial Gain used in calculation of the new xi 
Gamma=G; 
 
X=zeros(size(Y));                   % Initialize X (memory allocation) 
 
h=waitbar(0,'Wait...'); 
 
for i=1:numberOfFrames 
                                    %VAD and Noise Estimation START 
    if i<=NIS                        % If initial silence ignore VAD 
        SpeechFlag=0; 
        NoiseCounter=100; 
    else % Else Do VAD 
        [NoiseFlag, SpeechFlag, NoiseCounter, Dist]=vad(Y(:,i),N,NoiseCounter); 
    end                                                             % Magnitude Spectrum Distance VAD 
     
    if SpeechFlag==0                                                % If not Signal Update Noise Parameters 
        N=(NoiseLength*N+Y(:,i))/(NoiseLength+1);                   % Update and smooth noise mean 
        LambdaD=(NoiseLength*LambdaD+(Y(:,i).^2))./(1+NoiseLength); % Update and smooth noise variance 
    end 
    %VAD and Noise Estimation END 
     
    gammaNew=(Y(:,i).^2)./LambdaD;                       % A postiriori SNR 
    xi=alpha*(G.^2).*Gamma+(1-alpha).*max(gammaNew-1,0); % Decision Directed Method for A Priori SNR 
     
    Gamma=gammaNew; 
    nu=Gamma.*xi./(1+xi);                                 % A Function used in Calculation of Gain 
     
    G= (xi./(1+xi)).*exp(.5*expint(nu));                  % Log spectral MMSE [Ephraim 1985] 
    X(:,i)=G.*Y(:,i);                                     % Obtain the new Cleaned value 
     
    waitbar(i/numberOfFrames,h,num2str(fix(100*i/numberOfFrames))); 
end 
 
close(h); 
output = OverlapAdd2(X,YPhase,W,SP*W);                       % Overlap-add Synthesis of speech 
speech = filter(1,[1 -pre_emph],output);                     % Undo the effect of Pre-emphasis 

%% Feature Extraction 

% Feature extraction (feature vectors as columns)
[ MFCCs, FBEs, frames ] = ...
  mfcc( speech, fs, Tw, Ts, alpha, @hamming, [LF HF], M, C+1, L );
                
% Record all the features to the matrix
MFCCs = MFCCs(:);
FBEs = FBEs(:);
Test = [MFCCs; FBEs];
Test = Test(274:end);

% Load pre-trained fuzzy neural network
load Network;

% Feed the input to the network
output = sim(net,Test);
output = evalfis(output,fis);

% Round off the value to the nearest whole number
% microbes = 1,  Magnaporthe grisea
% microbes = 2,  Thanatephorus cucumeris
% microbes = 3,  Xanthomonas oryzae
microbes = round(output);
end