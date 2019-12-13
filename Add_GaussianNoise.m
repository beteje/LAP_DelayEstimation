function y_hat = Add_GaussianNoise(y, SNR)
%% adds White Gaussian Noise to input signal y based on the SNR value. Note
% SNR = 10*log10(sum(abs(y).^2)/sum(abs(y-y_hat).^2)

[M,N] = size(y);
e = randn(M,N);

y_hat = zeros(M,N);
for l = 1:M
    sigma = sqrt(norm(y(l,:)).^2*10^(-SNR/10)/norm(e(l,:)).^2);

    y_hat(l,:) = y(l,:) + e(l,:).*sigma;
end
