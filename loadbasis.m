function basis=loadbasis(K)
%% Generates the elements for the filter basis. 
% Based on the Gaussian filter basis defined in the ICASSP2015 and ICIP2015 papers.
%
% input:    K       - size of filter (2K+1)
%
% output:   basis   - Gaussian filter basis
%
% ICASSP2015: C. Gilliam and T. Blu, Local All-Pass Filters for Optical Flow Estimation
% ICIP2015: T. Blu, P. Moulin and C. Gilliam, Approximation Order of the LAP Optical Flow Algorithm

s = K/2-0.2;                                % Equation for sigma
g=@(x)exp(-x.*x/2/s^2);                     % Gaussian function

K0 = ceil(K);
k = (-K0:K0)';                              % Indices

basis = zeros(2*K0+1,2);
basis(:,1) = g(k); 
basis(:,1)=basis(:,1)/(sum(basis(:,1)));    % Normalize Gaussian
basis(:,2)=g(k).*k;                         % 1st derivative 
basis(:,2)=basis(:,2)/(sum(k.*basis(:,2))); % Normalize 1st derivative