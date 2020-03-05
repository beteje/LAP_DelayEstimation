function [theta,d_est,d_est_Kalman,d_est_Kalman_fus,MAE] = Delay_Est_Kalman(T,SNR,Fs,delta_e,vtype)
%% Estimates the delay between channels of data using both the LAP & Multiscale LAP
% Data is generated using a simple model which creates data with known spectrum
%
% inputs:   T                   - number of seconds of data to use
%           SNR                 - signal to noise ratio (dB)
%           Fs                  - sampling rate
%           delta_e             - distance between electrodes (mm)
%           vtype               - type of velocity 
%                                   1 = linear, 2 = sinusoidal, 3 = sigmoidal, otherwise constant
%
% outputs:  theta               - true delay signal
%           d_est               - LAP estimate of the delay
%           d_est_Kalman        - LAP+Kalman estimate of the delay
%           d_est_Kalman_fus    - LAP+Kalman fused estimate of the delay
%           MAE                 - mean absolute error of the delay estimates

K = [8,16];     % Half length of filter basis
q = 0.5;        % Process noise variance
r = 1;          % Measurement noise variance

dt = 1/Fs;      % Sampling interval
t = 0:dt:T;     % Time points
N_Sig = 2;      % Number of signals

% Generate Data
[x,theta] = Signal_Generation(length(t),N_Sig,SNR,Fs,vtype,4,0,0,1,delta_e,1,1,500); 

% Estimate the delays using LAP+Kalman
[d_est,d_est_Kalman,d_est_Kalman_fus] = LAP_Kalman(x,K,q,r,dt);

% Calculate MAEs
MAE(1:length(K)) = mean(abs(theta-d_est),2);
MAE(end+1:end+length(K)) = mean(abs(theta-d_est_Kalman),2);
MAE(end+1) = mean(abs(theta-d_est_Kalman_fus));