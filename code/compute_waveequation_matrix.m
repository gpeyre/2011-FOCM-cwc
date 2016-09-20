function [L,Delta] = compute_waveequation_matrix(c,derivative)

% compute_waveequation_matrix - compute the second order derivative matric
%
% [L,Delta] = compute_waveequation_matrix(c);
%
%   L=diag(c^2)*Delta is the matrix of c(x)^2*d^2/dx^2 (in 1D, but extends to 2D)
%
%   The domain is [0,1] or [0,1]^2 discretized.
%
%   Copyright (c) 2007 Laurent Demanet and Gabriel Peyre

if nargin<2
    derivative = 'fft';
end

ndims = 2;
if size(c,1)==1 || size(c,2)==1
    ndims = 1;
    c = c(:);
end
n = size(c,1);    

% compute 1D derivative matrix
Delta = compute_derivative_matrix(n, derivative);

if ndims==1
    L = diag(c.^2) * Delta;
else
    Delta = kron(eye(n),Delta)+kron(Delta,eye(n));
    L = diag(c(:).^2) * Delta;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = compute_derivative_matrix(n, derivative)

switch derivative
    case 'fd'
        h = 1/n;
        d2udx2 = zeros(n,1);
        d2udx2(1) = -2/h^2;
        d2udx2(2) = 1/h^2;
        d2udx2(end) = 1/h^2;
    case 'fft'
        L = 1;
        dx = L/n;

        % grid in x
        if mod(n,2)
            x = -(L/2-dx/2):dx:(L/2-dx/2);
        else
            x = -L/2:dx:(L/2-dx);
        end
        x = x';

        % grid in xi
        xi = 1:n;
        xifreq0 = floor(n/2+1);
        xi = xi - xifreq0;
        xi2 = xi.^2;
        xi = xi'; xi2 = xi2';

        u = zeros(n,1); u(1) = 1;

        Fdudx = (2*pi*1i*xi/L) .* fftshift(fft(ifftshift(u)));
        dudx = real(fftshift(ifft(ifftshift(Fdudx))));

        Fd2udx2 = (- 4*pi^2*xi2/L^2) .* fftshift(fft(ifftshift(u)));
        d2udx2 = real(fftshift(ifft(ifftshift(Fd2udx2))));
end

% d2udx2 = -2/h^2 * d2udx2/d2udx2(1);
[Y,X] = meshgrid(1:n,1:n);
D = d2udx2( mod(X-Y,n)+1 );
D = (D+D')/2; % ensure symetry
