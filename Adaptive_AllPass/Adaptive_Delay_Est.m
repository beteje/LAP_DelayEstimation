function [d,d_AP,d_Sun,d_ETDE,MSE,MAE] = Adaptive_Delay_Est(T,SNR,Fs,d1,step,mu)
%% Estimates the delay between channels of data using the NAAP, Sun Adaptive AllPass and ETDE
% Data is generated using a simple model which creates bandlimited data
%
% inputs:   T           - number of seconds of data to use
%           SNR         - signal to noise ratio (dB)
%           Fs          - sampling rate
%           d1          - the initial delay (in samples)
%           step        - the size of the largest step change
%           mu          - learning rates for each of the algorithms
%                           mu(1) - NAAP, mu(2) - Sun, mu(3) - ETDE
%
% outputs:  theta       - true delay signal
%           d_est       - LAP estimate of the delay
%           d_est_multi - multiscale LAP estimate of the delay
%           MAE         - mean absolute error of the delay estimates

%% Generate Data
t = 0:1/Fs:T;                           % Time points
N = length(t);                          % Number of samples
ds = [d1 d1+step/2 d1-step/2];          % Set the delay values
NumCon = length(ds);                    % Number of constant delays     
N1 = floor(N/NumCon);                   % Equal number of samples for each section
d = [];                                 % Generate the delay signal
for i=1:NumCon
    if i==NumCon
        d = [d ds(i)*ones(1,N-(i-1)*N1)];
    else
        d = [d ds(i)*ones(1,N1)];
    end
end
% Generate piecewise constant signal
[x,~] = Signal_Generation(N,2,SNR,Fs,5,ds);

%% Delay estimation
K = 7;              % Length of filters
con_start = 100;    % Set the number of samples before the constrained optimisation of ETDE

% Estimate delay with all-pass NLMS
[d_AP,e_AP,~] = Adaptive_AllPass(x.',K,mu(1));
MSE(1) = mean(abs(e_AP).^2);             % Calculate MSE
MAE(1) = mean(abs(d(1:end-K)-d_AP));     % Calculate MAE
% Estimate delay with Suns algorigthm
[d_Sun,e_Sun,~] = Adaptive_AllPass_Sun1999(x.',K,mu(2));
MSE(2) = mean(abs(e_Sun).^2);            % Calculate MSE
MAE(2) = mean(abs(d(1:end-K)-d_Sun));    % Calculate MAE
% Estimate delay with ETDE
[d_ETDE,e_ETDE] = ETDE(x.',2*K+1,mu(3),con_start);
MSE(3) = mean(abs(e_ETDE).^2);           % Calculate MSE
MAE(3) = mean(abs(d-d_ETDE(1:end-1)));   % Calculate MAE
