function [D,e] = ETDE(x_signal,K,mu,con_start)
%% Estimate delay between signals using explicit time delay estimation
% Based on "A New Algorithm fo Explicit Adaptation of Time Delay", So, Ching & Chan, 1994 IEEE TSP
%
% inputs:   x           - signals to estimate the delay between 
%           K           - length of the filter
%           mu          - step-size 
%           norm_trig   - trigger for the normalised version
%           w           - initial estimates of the filter weights
%
% outputs:  delay_est   - per sample estimate of the delay
%           e           - adaptive filter prediction error
%           w_holder    - estimates of the filter weights

if nargin == 2
    mu = 0.005;             % Set step-size
    con_start = 0;          % Start constrained optimisation immediately 
elseif nargin == 3
    con_start = 0;   
end

% Initialise parameters:
sample_size = length(x_signal);
P = floor((K-1)/2);

D = zeros(1,sample_size);
e = zeros(1,sample_size);
x = zeros(K, 1);
y = zeros(K, 1);
w = zeros(K, 1);
ind = -P:P;

% Adaptive algorithm
for i = 1:sample_size
    x(1) = x_signal(i,1);       % Select current samples
    y(1) = x_signal(i,2);
    
    % Initially use unconstrained optimisation to give integer estimation
    if i <= con_start
        e(i) = y(P+1) - w'*x;             
        w = w + mu*e(i)*x;
        [~,largest_position_index] = max(w);
        D(i+1) = largest_position_index - P - 1;
    else                        % Constrained optimisation
        sinc_td = sinc(ind - D(i));
        td = (cos(pi*(ind-D(i))) - sinc_td)./(ind - D(i));
        td(isnan(td)) = 0;
        e(i) = y(P+1) - sinc_td*x;
        D(i+1) = D(i) - mu*e(i)*td*x;
    end
    
    x(2:end) = x(1:end-1);
    y(2:end) = y(1:end-1);
end