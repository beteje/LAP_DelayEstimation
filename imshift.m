function J = imshift(I,u,kernel)
%% This function uses cubic OMOMS interpolation used in 'interp.m' with half point symmetric image extensions (can be changed).
% USAGE    : J=imshift(I,u);
% FUNCTION : Shifts the image I by real(u) along the lines and by imag(u)
% along the columns. 'u' is an arbitrary complex vector field. If 'u' is a 
% scalar, all the pixels are shifted by the same amount.
%
% DATE     : 23 November 2014
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org
%
% Related References:
% 1) T. Blu, P. Thevenaz, and M. Unser, "MOMS: Maximal-order interpolation 
%   of minimal support," IEEE Trans. Image Processing, vol. 10, no. 7, 
%   pp. 1069-1080, 2001.
% 2) T. Blu, P. Thevenaz, and M. Unser, "Linear interpolation revitalized,"
%   IEEE Trans. Image Processing, vol. 13, no. 5, pp. 710–719, 2004.

[M,N,P] = size(I);

if nargin == 2
    kernel='cubicOMOMS';
end

if length(u)>1              % Shift individual pixels by different amounts
    [M0,N0] = size(u);
    if M0==M && N0==N
        x = (1:M)'*ones(1,N);
        y = ones(M,1)*(1:N);
        J = zeros(M,N,P);
        for p=1:P
            J(:,:,p) = interp(x-real(u),y-imag(u),double(I(:,:,p)),kernel);
        end
        J = cast(J,'like',I);
    else
        error('Input image and flow dimensions do not match!')
    end
else                        % Shift all pixels by the same amount
    integer_shift = (abs(round(u)-u)<=1e-6);

    if integer_shift        % Shift all pixels by integer value
        u = round(u);
        xshift = abs(real(u));
        yshift = abs(imag(u));
        switch sign(real(u))
            case -1
                I1 = [I;I(end-(1:xshift)+1,:)];
            case 0
                I1 = I;
            case 1
                I1 = [I((xshift:-1:1),:);I];
        end
        switch sign(imag(u))
            case -1
                I1 = [I1 I1(:,end-(1:yshift)+1)];
            case 0
            case 1
                I1 = [I1(:,(yshift:-1:1)) I1];
        end

        J = I1((1:M)+max(0,-real(u)),(1:N)+max(0,-imag(u)));
    else                    % Shift all pixels by non-integer value
        x = (1:M)'*ones(1,N);
        y = ones(M,1)*(1:N);
        J = zeros(M,N,P);
        for p = 1:P
            J(:,:,p) = interp(x-real(u),y-imag(u),double(I(:,:,p)),kernel);
        end
        J = cast(J,'like',I);
    end
end