function coh = cohF_multi(x,segL)
%% Calculate coherence based on Welch averaging method
%
% inputs:   x1      - first signal
%           x2      - second signal
%           segL    - length of the data segment windows 
%                       creates 3 spectra of length segL/2 + 1 with 50% overlap
%
% outputs:  coh     - coherence

W = segL/2;         % Window length
L = segL/4;         % Step length
[M,N] = size(x);    % Original data length

% Create window
h = hann(W+1);

% Pad data
x = [fliplr(x(:,2:W+1)),x,fliplr(x(:,N-W:N-1))];

coh = zeros(W+1,N);
for n=W+1:N+W
    %Calculate windowed Fourier transforms of data
    X = zeros(W+1,3,M);
    for m=1:M
        for i=-1:1
            X(:,i+2,m) = fft(x(m,n+(i-1)*L:n+(i+1)*L).*h');
        end
    end
    
    Sx1x2 = zeros(W+1,1);
    Sx1x1 = zeros(W+1,1);
    Sx2x2 = zeros(W+1,1);
    % Calculate averaged spectra
    for m=1:M-1
        Sx1x2 = Sx1x2 + sum(X(:,:,m).*conj(X(:,:,m+1)),2);
        Sx1x1 = Sx1x1 + sum(X(:,:,m).*conj(X(:,:,m)),2);
        Sx2x2 = Sx2x2 + sum(X(:,:,m+1).*conj(X(:,:,m+1)),2);
    end
    
    % Calculate coherence of current data point
    coh(:,n-W) = Sx1x2./(sqrt(Sx1x1).*sqrt(Sx2x2));
end