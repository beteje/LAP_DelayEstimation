function [d_est,chosen_order,AIC] = MultiScale_LAP_Param(x1,x2,Level_Num,Min_Wind,Model_Order_array)
%% Implementation of a multi-scale framework for the LAP algorithm with parametric model fitting 
% The filter basis used in the LAP algorithm spans the derivatives of a Gaussian filter. Level_Num determines the number of scales used.
%
% inputs:   x1                  - input signal 1, Sig_Num (number of signals) by s1 array
%           x2                  - input signal 2, Sig_Num (number of signals) by s1 array
%           Level_Num           - number of scales used in the algorithm
%                                   e.g 5 scales = [32,16,8,4,2,1]
%           Min_Wind            - minimum size for the local window
%           Model_order_array   - Model orders to test
%                                   e.g. [1:3] = test models with orders 1, 2, and 3
%                                   and choose the best using AIC
%
% outputs:  d_est               - estimate of the time delay per sample
%           chosen_order        - the chosen model order
%           AIC                 - return the AIC used to select the model order

if nargin < 5
    % set maximum number of polynomial bases
    Num_coeff = 2;
    Model_Order_array = 1:Num_coeff;
else
    Num_coeff = max(Model_Order_array);
end

[Sig_Num, N] = size(x1);                        % Number of signals by length of signals
n = 1:N;                                        % Sample numbers

% time axis
t_index = linspace(-1,1,N);

% Chebyshev polynomial basis of first kind:
vec = cos(acos(t_index(:))*(0:(Num_coeff-1)));

% Initialise variables
d_est = zeros(1,N);                             % Output estimate of delay

if (2^(Level_Num+1)+1) > N
    warning(['Level Number of ', int2str(Level_Num), ' results in a filter larger than the size of the signal. Level number reduced']);
    Level_Num = floor(log2(N-1)-1);
end

% Generate array of filter sizes
LN = Level_Num+1;
amp_array = 2.^(linspace(Level_Num,1,floor(LN)));

% Estimate the delay
for l = 1:LN
    amp_size = amp_array(l);                    % Define scale factor  
    Basis_Set = loadbasis(amp_size);            % Load Filter Basis
    
    if l == 1 && sum(d_est) == 0
        % no warping on first iteration
        x2_shift = x2;
    else
        for index_sig = 1:Sig_Num
            % Using current delay estimate warp the electrodes in x2 closer to the electrodes in x1
            % Shift value needs to be imaginary to shift along x axis 
            x2_shift(index_sig,:) = imshift(x2(index_sig,:),-1i.*d_est); 
            
            % Sort out edge effects
            mask = logical((n + d_est)<=N & (n + d_est)>= 1);
            x2_shift(index_sig,mask == 0) = x1(index_sig,mask==0);
        end
    end
    
    wind = max([2*ceil(amp_size)+1,Min_Wind],[],2);         % Calculate local window size for LAP
    
    d_raw = LAP_1D(x1,x2_shift,Basis_Set,amp_size,wind).';  % Estimate delay using LAP
    
    % Remove nans using inpainting
    index1 = (round(wind)+1):(N-round(wind));
    d_clean = [d_raw(index1(1)).*ones(1,(index1(1)-1)), d_raw(index1), d_raw(index1(end)).*ones(1,(index1(1)-1))];
        
    % index of points to be used in the fitting:
    index1 = logical((round(wind)+1) <= (1:N) & (N-round(wind)) >= (1:N)).*~isnan(d_clean);
    
    % Update d_est:
    d_est = d_est + d_clean;    % Add estimated flow to the previous flow estimate
    
    AIC = NaN(length(Model_Order_array),1);
    d_fit = NaN(length(Model_Order_array),N);
    for poly_order1 = 1:length(Model_Order_array)
        poly_order = Model_Order_array(poly_order1);
        
        % Perform fitting:
        d_fit(poly_order1,:) = (vec(:,1:poly_order)*(vec(index1==1,1:poly_order)\d_est(index1==1).')).';
        
        % Calculate information theoretic critera:
        AIC(poly_order1) = 2*sum(abs(d_est - d_fit(poly_order1,:)).^2) + 2*N*poly_order/(N-1 - poly_order);           
    end
     
    [~,i2] = min(AIC);
    
    d_est = d_fit(i2,:);
    chosen_order = Model_Order_array(i2);
end