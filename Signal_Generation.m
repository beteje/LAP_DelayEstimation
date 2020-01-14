function [x,theta] = Signal_Generation(N,N_Sig,SNR,Fs,varargin)
%% Generate simple bandlimited multichannel data 
% Default is a badlimited signal with constant delay of 4 samples
%
% inputs:   N           - number of data points to generate
%           N_Sig       - number of signals to generate
%           SNR         - signal to noise ratio
%                           either a single value or an array of values
%           Fs          - sampling frequency
%
% optional inputs:
%           type        - type of delay/velocity
%                           1 = linear (delay), 2 = sinusoidal, 3 = sigmoidal, 
%                           4 = smoothly varying bandlimited, 
%                           otherwise constant (default)
%           d1          - the constant part of the delay
%                           default = 4
%           alpha       - Defines the bandwidth of the signal as B = Fs*alpha
%                           default = 0.5
%           pm          - Number of coeffficients for interpolation (2pm + 1)
%                           0 = use OMOMS (default)
%           velocity    - define the delay in terms of velocity for sinusiodal, sigmoidal or constant options
%                           1 = Yes, otherwise no (default)
%           delta_e     - electrode spacing (used to caluclate velocity)
%                           default = 0.005
%           psd         - define a specified power spectral density
%                           1 = Yes, otherwise no (default)
%           fltSig      - filter the noisy signals 
%                           1 = Yes, otherwise no (default)
%           lpf         - frequency for the low-pass filter
%                           default = 500Hz
%
% outputs:  x       - synthetic signals (N_Sig by N)

% Check number of inputs
narginchk(4,13);

% Set default parameters
defaultType = 0;
defaultD1 = 4;
defaultAlpha = 0.5;
defaultPm = 0;
defaultVelocity = 0;
defaultDelta_e = 0.005;
defaultPSD = 0;
defaultFltSig = 0;
defaultLPF = 500;

% Create parser
p = inputParser;
addOptional(p,'type',defaultType);
addOptional(p,'d1',defaultD1);
addOptional(p,'alpha',defaultAlpha);
addOptional(p,'pm',defaultPm);
addOptional(p,'velocity',defaultVelocity);
addOptional(p,'delta_e',defaultDelta_e);
addOptional(p,'psd',defaultPSD);
addOptional(p,'fltSig',defaultFltSig);
addOptional(p,'lpf',defaultLPF);
% Parse variable input arguments
parse(p,varargin{:});

x = zeros(N_Sig,N);                     % Initialize signals   
f = linspace(-Fs/2,Fs/2,N);             % Define frequency points

% % lpf = 500;                          % Frequency to low pass filter noisy data at

% Set type of delay/velocity
if p.Results.type==1            % linear 
    d = p.Results.d1 + 2*linspace(0,1,N);
elseif p.Results.type==2        % sinusoidal
    d = p.Results.d1 + 2*sin(2*pi*0.2*(1:(N))/Fs);
elseif p.Results.type==3        % sigmoidal
    d = p.Results.d1 + 1./(1 + exp(-8/5.*((0:(N-1))./Fs - (N-1)/Fs./4)));
elseif p.Results.type==4        % smoothly varying bandlimited   
    d = randn(1,N);
    freq = linspace(-Fs/2,Fs/2,N);
    Filter = ifftshift(double(logical(abs(freq)<=10)));
    d = real(ifft(fft(d).*Filter));
    d = d./(max(abs(d))).*6;    
else                            % constant
    d = p.Results.d1.*ones(1,N);
end

% If calculating velocity convert d into delay
% (cannot use velocity option for linear or bandlimited delay)
if p.Results.velocity == 1 && p.Results.type~=4
    theta = Fs*p.Results.delta_e./d; 	% d is velocity, theta is delay
else
    theta = d;                          % d is delay
end

% Select whether to use the defined power spectral density or a bandlimited signal
if p.Results.psd ==1   
    fl = 60;                            % Low signal frequency      
    fh = 120;                           % High signal frequency
    P = fftshift((1*fh^4*f.^2)./((f.^2+fl^2).*(f.^2+fh^2).^2)); % Desired PSD  
else
    B = p.Results.alpha.*Fs./2;  % Define bandwidth of the signals
    bwf = zeros(1,N);
    bwf(:,find(f==-B):find(f==B)) = 1;  % Create brick wall filter with desired cutoff
    P = ifftshift(bwf);
end

X1 = fft(randn(1,N));                   % Frequency response of noise
X1 = X1.*P;                             % Frequency response with desired PSD
x(1,:) = real(ifft(X1));                % First simulated channel

i = -p.Results.pm:p.Results.pm;         % Indices
% Iterate through to calculate subsequent channels
for j=2:N_Sig
    if p.Results.pm==0
        % Use cubic OMOMS for interpolation
        x(j,:) = imshift(x(j-1,:),1i.*theta);
    else
        % Use sinc funciton for interpolation
        tmp = [fliplr(x(j-1,2:p.Results.pm+1)),x(j-1,:),fliplr(x(j-1,N-p.Results.pm:N-1))];
        for n=1+p.Results.pm:N+p.Results.pm
            w = sinc(i-theta(n-p.Results.pm));        % Create filter coefficients for current time point
            x(j,n-p.Results.pm) = w*tmp(n-i).';       % Create new signal as delayed version of previous
        end
    end
end

x = Add_GaussianNoise(x,SNR);           % Add noise to channels
if p.Results.fltSig == 1
    % Create brick wall filter with desired cutoff
    bwf = zeros(N_Sig,N);
    bwf(:,find(f==-p.Results.lpf):find(f==p.Results.lpf)) = 1;  
    X = fft(x,[],2);                    % Take frequency response of noisy channels
    X = X.*fftshift(bwf);            	% Filter the channels
    x = real(ifft(X,[],2));          	% Convert filtered channels back to time
end
