function [theta,d_est,d_est_multi,theta_cohF,MAE] = Delay_Est(N_Sig,T,SNR,Fs,delta_e,cv_type)
%% Estimates the delay between channels of data using both the LAP & Multiscale LAP
% Data is generated using a simple model which creates data with known spectrum
% CohF function is also included as a comparison
%
% inputs:   N_Sig       - number of signals to generate
%           T           - number of seconds of data to use
%           SNR         - signal to noise ratio (dB)
%           Fs          - sampling rate
%           delta_e     - distance between electrodes (mm)
%           cv_type     - type of conduction velocity 
%                           1 = linear, 2 = sinusoidal, 3 = sigmoidal, otherwise constant
%
% outputs:  theta       - true delay signal
%           d_est       - LAP estimate of the delay
%           d_est_multi - multiscale LAP estimate of the delay
%           MAE         - mean absolute error of the delay estimates

%% Generate Data
t = 0:1/Fs:T;	% Time points
% Generate sample channels
[x,theta] = Simple_EMG_Model(length(t),N_Sig,SNR,Fs,delta_e,cv_type); 

%% LAP delay estimation
% Half length of filter basis
K = 11;         % Minimum value to estimate a delay tau is K = 2*tau;
% Window length
W = 513;        % Minimum value is W = 2*K + 1 (larger values improve performance under noise)                                

% Generate filter basis
basis = loadbasis(K);

% Estimate per sample delay using LAP
d_est = LAP_1D(x(1:N_Sig-1,:),x(2:N_Sig,:),basis,K,W);
MAE(1) = nansum(abs(theta-d_est'))/length(t);   % Calculate MAE

%% Multiscale LAP delay estimation
N_Scales = 3;       % Number of scales to estimate across
Min_Wind = 513;     % Minimum window size

% Estimate per sample delay using multiscale Lap
d_est_multi = MultiScale_LAP(x(1:N_Sig-1,:),x(2:N_Sig,:),N_Scales,Min_Wind);
MAE(2) = sum(abs(theta-d_est_multi))/length(t); % Calculate MAE

%% Coherence delay estimation - spectrum averaging
cohF_est = cohF_multi(x,Fs/2);                  % Calculate spectrum averaging coherence
f_incF = 4;                                     % Frequency increments for cohF

f_min = 16;                                     % Minimum frequency of interest
f_max = 200;                                    % Maximum frequency of interest

f_cohF = (f_min:f_incF:f_max).';                % Frequencies for cohF
indF = (f_min/f_incF+1:f_max/f_incF+1);         % Indices for cohF frequencies

angle_cohF = angle(cohF_est(indF,:));           % Phase angles from cohF
theta_cohF = f_cohF\angle_cohF;                 % Linear fit of delay at each time point
theta_cohF = theta_cohF*Fs/(2*pi);           	% Adjust delay by Fs/2pi

MAE(3) = sum(abs(theta-theta_cohF))/length(t);	% Calculate mean absolute error