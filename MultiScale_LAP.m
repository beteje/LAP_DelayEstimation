function d_est = MultiScale_LAP(x1,x2,Level_Num,Min_Wind)
%% Implementation of a multi-scale framework for the LAP algorithm. 
% The filter basis used in the LAP algorithm spans the derivatives of a Gaussian filter. Level_Num determines the number of scales used.
%
% inputs:   x1          - input signal 1, Sig_Num (number of signals) by s1 array
%           x2          - input signal 2, Sig_Num (number of signals) by s1 array
%           Level_Num   - number of scales used in the algorithm
%                           e.g 5 scales = [32,16,8,4,2,1]
%           Min_Wind    - minimum size for the local window
% outputs:  d_est       - estimate of the time delay per sample

[Sig_Num, N] = size(x1);                        % Number of signals by length of signals
n = 1:N;                                        % Sample numbers

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
    d_raw = [d_raw(index1(1)).*ones(1,(index1(1)-1)), d_raw(index1), d_raw(index1(end)).*ones(1,(index1(1)-1))];
    d_clean = clean_delay(d_raw, not(isnan(d_raw)));
    
    d_clean = Cleaning_Procedure_Mean(d_clean, amp_size);   % Smooth flow estimation
    
    d_est = d_est + d_clean;    % Add estimated flow to the previous flow estimate
end
end

function d_out = Cleaning_Procedure_Mean(d_in, amp_size)
%% Cleans the delay estimate. 
% The cleaning process comprises smoothing using a Gaussian filter

R1 = round(2*amp_size);             % Sigma value for Gaussian filter

% Obtain filters
k1 = (-R1:R1);
Gauss_Filter = exp(-k1.*k1/2/R1^2);
Gauss_Filter = Gauss_Filter./norm(Gauss_Filter(:),1);

d_out = imfilter(d_in, Gauss_Filter,'symmetric');   % Smooth input
end