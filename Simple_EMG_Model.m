function [x,theta] = Simple_EMG_Model(N,N_Sig,SNR,Fs,delta_e,type)
%% Generate simple multichannel EMG like data with known spectrum 
% Assumes spatially separated electrodes with the delay based on the conduction velocity
%
% inputs:   N       - number of data points to generate
%           N_Sig   - number of signals to generate
%           Fs      - sampling frequency
%           delta_e - electrode spacing
%           type    - type of MFCV
%                       1 = linear, 2 = sinusoidal, 3 = sigmoidal, otherwise constant
%
% outputs:  x       - simulated EMG channels (N_Sig by N)
%                       each channel after 1st is a delayed version of precdeding

x = zeros(N_Sig,N);                 % Initialize signals            

fl = 60;                            % Low signal frequency      
fh = 120;                           % High signal frequency
pm = 0;                             % Number coefficients for interpolation (2pm+1)
                                        % if pm=0 use OMOMS
lpf = 500;                          % Frequency to low pass filter noisy data at

f = linspace(-Fs/2,Fs/2,N);         % Frequency points
bwf = zeros(N_Sig,N);
bwf(:,find(f==-lpf):find(f==lpf)) = 1;% Create brick wall filter with desired cutoff

% Set type of muscle fibre conduction velocity
if type==1          % linear 
    d = 2/delta_e + 1.5/delta_e*linspace(0,1,N);
    MFCV = 1./d;
elseif type==2      % sinusoidal
    MFCV = 4 + 2*sin(2*pi*0.2*(1:(N))/Fs);
elseif type==3      % sigmoidal
    MFCV = 2 + 1./(1 + exp(-8/5.*((0:(N-1))./Fs - (N-1)/Fs./4)));
else                % constant
    MFCV = 4;
end

theta = delta_e./MFCV;           % Time delay

P = fftshift((1*fh^4*f.^2)./((f.^2+fl^2).*(f.^2+fh^2).^2)); % Desired PSD 

X1 = fft(randn(1,N));               % Frequency response of noise
X1 = X1.*P;                         % Frequency response with desired PSD
x(1,:) = real(ifft(X1));            % First simulated channel

i = -pm:pm;                         % Indices
% Iterate through to calculate subsequent channels
for j=2:N_Sig
    if pm==0
        % Use cubic OMOMS for interpolation
        x(j,:) = imshift(x(j-1,:),1i.*theta);
    else
        % Use sinc funciton for interpolation
        tmp = [fliplr(x(j-1,2:pm+1)),x(j-1,:),fliplr(x(j-1,N-pm:N-1))];
        for n=1+pm:N+pm
            w = sinc(i-theta(n-pm));        % Create filter coefficients for current time point
            x(j,n-pm) = w*tmp(n-i).';       % Create new signal as delayed version of previous
        end
    end
end

x = Add_GaussianNoise(x,SNR);           % Add noise to channels
X = fft(x,[],2);                        % Take frequency response of noisy channels
X = X.*fftshift(bwf);                   % Filter the channels
x = real(ifft(X,[],2));                 % Convert filtered channels back to time
