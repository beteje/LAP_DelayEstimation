function [d_est] = LAP_1D(x1,x2,basis,K1,W,Raw)
%% Implementation of the Local All-Pass Filter (LAP) algorithm for delay estimation in 1D signals 
% Given input signals I1 and I2, and basis (the set of filters), the function calculates the delay, possibly time-varying, between the two signals.
%
% inputs:   x1      - input signal 1, Sig_Num (number of signals) by s1 array
%           x2      - input signal 2, Sig_Num (number of signals) by s1 array
%           basis   - elements of the filter basis. M by N array 
%                       N = number of basis elements 
%                       M = R^2 where the filter is R by R in size.
%           K1      - size of filter basis (Default = calculated from filter basis)
%           W       - size of the local window used in the algorithm.
%                       (Default = 2*ceil(K1) + 1)
%           Raw     - control variable to decide whether to keep raw estimate (1) or not (0)
%
% output:   d_est   - estimate of the time delay per sample

[Sig_Num, s1] = size(x1);   % Number of signals by length of signals
[M,N] = size(basis);        % Number of basis elements and square of filter size 

% Set default parameters
if nargin < 4
    K1  = (M-1)/2;
    W   = 2.*ceil(K1)+1;
    Raw = 0;
elseif nargin < 5
    W   = 2.*ceil(K1)+1;
    Raw = 0;
elseif nargin < 6
    Raw = 0;
end

% Check if basis filters are centered
K0  = 2*ceil(K1)+1;
K   = (K0-1)/2;
if round(K)~=K
    error('Basis filters are not centered.')
end

% Solving for input basis filters 
A = zeros(N-1,N-1,s1);
b = zeros(N-1,s1);
for sig_index = 1:Sig_Num
    II=[];
    % Filter each signal with the filter basis
    for n=1:N
        B1  = basis(:,n);
        B2  = B1(end:-1:1,1);
        J   = imfilter(x2(sig_index,:),B2.','symmetric')-imfilter(x1(sig_index,:),B1.','symmetric');
        II  = [II J(:)];
    end
    
    % Generate matrices equating to a linear systems of equations at each sample
    J = II(:,1);
    for k=1:N-1
        for l=k:N-1
            A(k,l,:) = A(k,l,:) + reshape(average(II(:,k+1).*II(:,l+1),W), 1,1,length(II(:,k+1)));
            A(l,k,:) = A(k,l,:);
        end
        b(k,:) = b(k,:) - shiftdim(average(II(:,k+1).*J,W),1);
    end
end


% Perform Gauss elimination on all samples in parallel
coeffs = zeros(s1,N-1);
for k=1:(N-1)
    for l=(k+1):(N-1)
        c = A(l,k,:)./A(k,k,:);
        for m=(k+1):(N-1)
            A(l,m,:) = A(l,m,:)-c.*A(k,m,:);
        end
        A(l,k,:) = 0;
        b(l,:) = b(l,:)-shiftdim(c,1).*b(k,:);
    end
end
for k=(N-1):-1:1
    coeffs(:,k) = shiftdim(b(k,:));
    for m=(k+1):(N-1)
        coeffs(:,k) = coeffs(:,k)-shiftdim(A(k,m,:)).*coeffs(:,m);
    end
    coeffs(:,k) = coeffs(:,k)./shiftdim(A(k,k,:));
end

if Raw ~= 1
    % Exclusion of the boundaries
    p = subindex(x1(1,:),K);
    coeffs(p,:) = NaN;
end

% First coefficient is always 1
coeffs = [ones(s1,1) coeffs];

% Determine delay estimates from the local all-pass filters
k0  = (-K:K)';
u1  = zeros(s1,1);
u11 = zeros(s1,1);
for n=1:N
    u1(:)   = u1(:)-(basis(:,n)'*k0(:))*coeffs(:,n);
    u11(:)  = u11(:)+sum(basis(:,n))*coeffs(:,n);
end
% per sample delay:
d_est = 2*(u1./(u11));

if Raw ~= 1
    % if delay is larger than the filter size then the estimate is erroneous
    % and set to nan
    d_est(abs(d_est)>K) = NaN;
end

%% Embedded functions 
function J=average(I,W)
%% Performs averaging over local region
J = imfilter(I.',ones(1,W),'symmetric')/(W);
J = J(:);
return

function p=subindex(I,K)
%% Determines boundary samples
N = length(I);
n = (1:N);
p = find(n<(2*K+1) | n> (N-2*K));
return