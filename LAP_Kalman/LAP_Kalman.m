function [d_est,d_est_Kalman,d_est_Kalman_fus] = LAP_Kalman(x,K,q,r,dt)
%% Implementation of the LAP algorithm with Kalman filter and fusion of multiple LAP+Kalman filters
%
% inputs:   x                   - input signals
%           K                   - size of filter basis for the LAP
%           q                   - process noise variance
%           r                   - measurement noise variance
%           dt                  - sampling interval
%
% outputs:  d_est               - estimate of the time delay per sample from the LAP for each K
%           d_est_Kalman        - estimate of the time delay from the LAP+Kalman for each K
%           d_est_Kalman_fus    - fused estimate of the time delay from the LAP+Kalman

A = [1 dt dt^2/2; 0 1 dt; 0 0 1];   % State transition
G = [dt^2/2; dt; 1];                % Process noise model
Q = G*q*G';                         % Process noise
C = [1 0 0];                        % Observation model

d_est = zeros(length(K),length(x));
d_est_Kalman = zeros(length(K),length(x));
for k = 1:length(K)
	% LAP delay estimation
    W = K(k);                 
    % Generate filter basis
    basis = loadbasis(K(k));
    d_est(k,:) = LAP_1D(x(1,:),x(2,:),basis,K(k),W,1)';
    
    % Individu
    x_tmp = indKalman(d_est(k,:),A,C,Q,r);
    d_est_Kalman(k,:) = x_tmp(1,:);
end

x_tmp = fusKalman(d_est,A,Q,r);
d_est_Kalman_fus = x_tmp(1,:);

end

%% Individual Kalman Filter
function x_est = indKalman(y,A,C,Q,R)
    I = eye(3);
    P = I;                                          % initial estimate error covariance
    x_est = zeros(3,length(y));                     % state estimate
    x_est(1,1) = y(1);
    At = A';
    for n = 2:length(y)
        x_pred = A*x_est(:,n-1);                    % predicted state estimate
        P_pred = A*P*At + Q;                        % predicted error covariance
        K = P_pred(:,1)/(R + P_pred(1,1));          % Kalman gain
    
        x_est(:,n) = x_pred + K*(y(n) - x_pred(1));	% updated state estimate
        P = (I-K*C)*P_pred;                         % updated estimated covariance
    end
end

%% Fused Kalman Filter
function x_est = fusKalman(y,A,Q,R)
    NStates = length(A);
    C_fus = repmat([1 zeros(1,NStates-1)],size(y,1),1); % fused observation model 
    P = eye(NStates);                                   % initial estimate error covariance
    x_est = zeros(NStates,length(y));                   % state estimate fusion 
    R_fus = diag(repmat(R,1,size(y,1)));                % fused covariance matrix 
    for n = 2:length(y) 
        % measurement fusion
        x_pred = A*x_est(:,n-1);                        % predicted state estimate
        P_pred = A*P*A' + Q;                            % predicted error covariance
        
        z_fus = y(:,n);                                 % measurement data
        K = P_pred*C_fus'/(R_fus + C_fus*P_pred*C_fus');	% Kalman gain
        
        x_est(:,n) = x_pred + K*(z_fus - C_fus*x_pred);	% updated state estimate
        P = (eye(NStates)-K*C_fus)*P_pred;              % updated estimated covariance
    end
end