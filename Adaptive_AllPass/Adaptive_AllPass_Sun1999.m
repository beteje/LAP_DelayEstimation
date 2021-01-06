function [delay_est,e,w_holder] = Adaptive_AllPass_Sun1999(x,K,mu)
%% Estimate delay between signals using an adaptive all-pass filter
% Based on "Adaptive Time Delay Estimation with AllPass Constraints", Sun & Douglas, Asilomar 1999
%
% inputs:   x           - signals to estimate the delay between 
%           K           - length of the filter
%           mu          - step-size 
%
% outputs:  delay_est   - per sample estimate of the delay
%           e           - adaptive filter prediction error
%           w_holder    - estimates of the filter weights

if nargin == 2
    mu = 0.05;      % Step-size
end

% initialise parameters:
x = [zeros(K,2); x];
w = zeros(K+1,1);
w_holder = zeros(K+1,length(x));
delay_est = zeros(1,length(x));
e = zeros(1,length(x));

index = toeplitz(K:(2*K),K:-1:0);

for n = (2*K+1):(length(x))
    % Samples n-K to n-2K of the first signal (x in paper)
    x1 = x((n-K):-1:(n-2*K),1); 
    
    % Matrix of past of second signal (z in paper)
    Z = reshape(x((n-index),2),K+1,K+1);
    
    % Filtered signal y at n-K
    y = w'*x1;
    
    % z at n-K, scalar:
    z = x(n-K,2);
    
    % Error term
    e(n) = z - y;
    
    % Calculate signal u in the paper:
    u = sum(Z.*(w_holder(:,(n-1)-(0:K)).'),2);
    
    % Update of w:
    w = w + mu.*(z*x1 - y*u);   
    w_holder(:,n) = w;
   
    % Estimate the delay
    w1 = [0; w; 0];
    [~,i2] = max(w1);
    i2 = max([2,i2]);
    if w1(i2-1) < w1(i2+1)
        delay_est(n) = (i2-2) + w1(i2+1)/(w1(i2+1)+w1(i2));
    else
        delay_est(n) = (i2-2) - w1(i2-1)/(w1(i2-1)+w1(i2));
    end
    
end

delay_est = delay_est((2*K+1):end);
e = e((2*K+1):end);
w_holder = w_holder(:,(2*K+1):end);
