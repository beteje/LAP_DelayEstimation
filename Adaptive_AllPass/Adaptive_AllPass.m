function [delay_est,e,w_holder] = Adaptive_AllPass(x,K,mu,norm_trig,w)
%% Estimate delay between signals using an adaptive all-pass filter
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
    mu = 0.05;              % Set step-size
    norm_trig = 1;          % Select normalised
    w = zeros(K,1);         % Initialise weights to 0
elseif nargin == 3
    norm_trig = 1;
    w = zeros(K,1);
elseif nargin == 4
    w = zeros(K,1);    
end

% Initialise parameters:
w_holder = w;
delay_est = zeros(1,length(x));
e = zeros(1,length(x));

holder = 1:K;
x1 = zeros(K,1);

for n = 1:(length(x)-K)
    if n<(K+2)    
        index = n-holder;
        index1 = sort(index(index>0));
        x1(index1) = x(index(index1),1);                    % Select the K past samples first data set
    else
        x1 = x(n-1:-1:n-K,1);                               % Select the K past samples first data set
    end
    x2 = x(n+1:n+K,2);                                      % Select the K future samples of the second data set
    
    e(n) = x(n,2) - x(n,1) + w'*(x2 - x1);                  % Calculate the error including the current data points
    
    if norm_trig == 0
        w = w - mu*e(n)*(x2-x1);                            % Update the filter weights
    else
        w = w - mu*e(n)*(x2-x1)./(norm((x2-x1))^2+0.00001);
    end
    
    w_holder = [w_holder, w];                               % Store the filter weights
          
    delay_est(n)= 2*(((1:K)*w)/(1+sum(w)));                 % Calculate the delay estimate
end

delay_est = delay_est(1:n);
e = e(1:n);