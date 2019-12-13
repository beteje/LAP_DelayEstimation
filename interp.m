function [J,offset]=interp(x,y,I,inttype)
%% Image interpolation function
% USAGE    : [J,offset]=interp(x,y,I,[inttype])
% FUNCTION : Interpolates the image I at the (in general, noninteger) 
% positions x and y (usual Matlab convention for images). x and y are 2D 
% matrices with same size. 

% DEFAULT  : The optional argument inttype can take the values
%               * 'nearestneighbor'
%               * 'bilinear' (default)
%               * 'keys'
%               * 'cubicspline'
%               * 'cubicOMOMS'
%               * 'shiftedlinear'
% Assumes images in double precision.
%
% DATE     : 23 November 2014
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org
%
% Related References:
% 1) T. Blu, P. Thevenaz, and M. Unser, "MOMS: Maximal-order interpolation 
%   of minimal support," IEEE Trans. Image Processing, vol. 10, no. 7, 
%   pp. 1069-1080, 2001.
% 2) T. Blu, P. Thevenaz, and M. Unser, "Linear interpolation revitalized,"
%   IEEE Trans. Image Processing, vol. 13, no. 5, pp. 710�719, 2004.
%

if nargin<4
    inttype='bilinear';
end

[a,b]=size(I);

% Determination of the interpolation function
switch inttype
    
    case 'nearestneighbor'
        L1  = -0.5;    % Support of the interpolation kernel
        L2  = +0.5;
        phi = @nearest;
        
    case 'bilinear'
        L1  = -1;      
        L2  = +1;
        phi = @linspline;
        
    case 'keys'
        L1  = -2;      
        L2  = +2;
        phi = @keys;
        
    case 'cubicspline'
        L1  = -2;      
        L2  = +2;
        phi = @cubicspline;
        
    case 'cubicOMOMS'
        L1  = -2;      
        L2  = +2;
        phi = @cubicOMOMS;
        
    case 'shiftedlinear'
        tau = 1/2*(1-sqrt(3)/3);
        L1  = floor(-1+tau);      
        L2  = ceil(1+tau);
        phi = @(x)linspline(x-tau);        
end

% Minimum and maximum row index needed in the interpolation formula
k0 = floor(min(x(:))-L2+1);
k1 = floor(max(x(:))-L1);
l0 = floor(min(y(:))-L2+1);
l1 = floor(max(y(:))-L1);

offset = [1-k0 1-l0];

% Smallest box enclosing the image and the (x,y) positions
kk0 = min(k0,1);
kk1 = max(k1,a);
ll0 = min(l0,1);
ll1 = max(l1,b);

% Indices used in the interpolation formula
k = floor(x-L2+1);
l = floor(y-L2+1);

% Image extension to fill the unknown pixels
exttype = 'symh';                               % options are: 
                                                % 'zpd' (zero-padding), 
                                                % 'symh' (half-point symmetry), 
                                                % 'symw' (whole-point symmetry), 
                                                % 'ppd' (periodization)
I0      = ext(I,exttype,[1-kk0 kk1-a 1-ll0 ll1-b]);
I0      = I0(1-kk0+(k0:k1),1-ll0+(l0:l1));
[a0,~] = size(I0);

% Prefiltering when needed
switch inttype
    
    case 'cubicspline'
        J = symfilter(2/3,1/6,I);
        J = symfilter(2/3,1/6,J.').';
        I0 = ext(J,exttype,[1-kk0 kk1-a 1-ll0 ll1-b]);
        I0 = I0(1-kk0+(k0:k1),1-ll0+(l0:l1));

    case 'cubicOMOMS'        
        J = symfilter(13/21,4/21,I);
        J = symfilter(13/21,4/21,J.').';
        I0 = ext(J,exttype,[1-kk0 kk1-a 1-ll0 ll1-b]);
        I0 = I0(1-kk0+(k0:k1),1-ll0+(l0:l1));

    case 'shiftedlinear'
        z0 = tau/(1-tau);
        % along columns first
        I0 = 1/(1-tau)*filtering(1,[1 z0],I0,'causal');
        % then lines
        I0 = I0.';
        I0 = 1/(1-tau)*filtering(1,[1 z0],I0,'causal');
        I0 = I0.';
end

% Kernel-based interpolation formula
J = zeros(size(x));
for dk = 0:(L2-L1-1)
    for dl = 0:(L2-L1-1)
        ind = k + dk + offset(1) + a0*(l+dl+offset(2)-1);       % matrices are ordered along columns in Matlab
        J = J + phi(x-k-dk).*phi(y-l-dl).*I0(ind);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Embedded functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = ext(I,exttype,extsize)
%% Image Extension
[a,b]   = size(I);
newa    = a + extsize(1) + extsize(2);
newb    = b + extsize(3) + extsize(4);

if extsize(1)>extsize(2)
    J = wextend(2,exttype,I,extsize(1),'bn');
    J = J(1:newa,:);
elseif extsize(2)>extsize(1) || (extsize(2)==extsize(1) && extsize(1)~=0)
    J = wextend(2,exttype,I,extsize(2),'bn');
    J = J(end+(1:newa)-newa,:);
else
    J = I;
end

if extsize(3)>extsize(4)
    J = wextend(2,exttype,J,extsize(3),'nb');
    J = J(:,1:newb);
elseif extsize(4)>extsize(3) || (extsize(4)==extsize(3) && extsize(3)~=0)
    J = wextend(2,exttype,J,extsize(4),'nb');
    J = J(:,end+(1:newb)-newb);
else
    J = J;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = nearest(x)
%% Nearest Neighbor Kernel
u = (x<0.5&x>=-0.5);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = linspline(x)
%% Linear Spline Kernel
u = 1 - abs(x);
u = u.*(u>=0);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = keys(x)
%% Keys Kernel
a   = -1/2;
x   = abs(x);
x2  = x.*x;
x3  = x.*x2;
u   =((a+2)*x3-(a+3)*x2+1).*(x<=1)+...
        (a*x3-5*a*x2+8*a*x-4*a).*(x>1&x<=2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = cubicspline(x)
%% Cubic Spline Kernel
xi  = 2 - abs(x);
xi2 = xi.*xi;
xi3 = xi.*xi2;
u   = (2/3-2*xi+2*xi2-1/2*xi3).*(xi>1&xi<=2)+...
        (1/6*xi3).*(xi>0&xi<=1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = cubicOMOMS(x)
%% Cubic OMOMS Kernel
x1=1-abs(x);
x13=x1.*x1.*x1;
x2=2-abs(x);
x23=x2.*x2.*x2;
u=(1/6*x23-2/3*x13+1/42*x2-2/21*x1).*(x2>1&x2<=2)+...
    (1/6*x23+1/42*x2).*(x2>0&x2<=1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function J = filtering(b,a,I,type)
%% Pre-filtering
switch type
    case 'causal'
        J=bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:))),I(1,:)*sum(b)/sum(a));
    case 'anticausal'
        I = I(end:-1:1,:);
        J = bsxfun(@plus,filter(b,a,bsxfun(@minus,I,I(1,:))),I(1,:)*sum(b)/sum(a));
        J = J(end:-1:1,:);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = symfilter(a,b,x)
%% All-pole IIR filter
% USAGE    : y=symfilter(a,b,x)
% FUNCTION : Implements the IIR filtering with an all-pole filter with 
% z-transform of the form 1/(a+b(z+1/z)) where a and b are such that the 
% denominator has no roots on the unit circle. This is equivalent to having
% a/b complex-valued (not real) or |a/b|>2.
%
% DATE     : 19 December 2014
% AUTHOR   : Thierry Blu, mailto:thierry.blu@m4x.org

[N,~]=size(x);

% Find the root z0 that is inside the unit circle
z0 = roots([b a b]);
if abs(z0(1))<1
    z0 = z0(1);
else
    z0 = z0(2);
end
if abs(z0)>=1
    error('The filter has poles on the unit circle!')
end
A=1/b/(z0-1/z0);

% One-pole IIR filtering of a symmetrized version of x 
z0n = filter(1,[1 -z0],[1;zeros(2*N,1)]);
y   = filter(1,[1 -z0],[x;x(end:-1:1,:);x(1,:)]);
a   = (y(end,:) - y(1,:))./(1 - z0n(end));
y   = y + z0n*a;
y   = A*(y(1:N,:) + y(2*N+1-(1:N),:) - x);
return
