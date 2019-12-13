function d0=clean_delay(d,mask)
%% Inpainting procedure to remove NaN values from delay
% inputs:   d       - input delay signal
%           mask    - a signal with same size as d equal to
%                       * 0 for the samples of d to be inpainted
%                       * 1 for the samples of d to be kept as is
% 
% outputs:  d0      - cleaned delay signal

mask = double(mask);
d(mask~=1) = 0;

% Diffusion kernel:
h = [1,0,1];

d0 = d;
encore = 1;
mask0 = mask;
while encore
    dn = imfilter(mask0,h,'symmetric');
    f = imfilter(mask0.*d0,h,'symmetric');
    n = find(mask==0&dn~=0);
    d0(n) = f(n)./dn(n);
    mask0 = double(dn~=0);
    encore = double(min(dn(:))==0);
end
