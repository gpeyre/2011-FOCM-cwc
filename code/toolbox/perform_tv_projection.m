function [x,err_tv,err_l2] = perform_tv_projection(x0,tau,options)

% perform_tv_correction - perform correction of the image to that it minimizes the TV norm.
%
%   [x,err_tv,err_l2] = perform_tv_projection(x0,tau,options);
%
%   Perform the following projection on the TV ball of radius tau:
%       min_x |x-x0|_2   s.t.   |x|_TV <= tau
%
%   The method use a subgradient projection, as explained in:
%       P.L. Combettes, J.C. Pesquet
%       "Image restoration subject to a total variation constraint"
%       IEEE Transactions on Image Processing 13 nç9 (2004), p. 1213-1222
%   The equations for the projection on the intersection of 2 half space is in
%       P.L. Combettes,
%       "A Block-Iterative Surrogate Constraint Splitting Method for Quadratic Signal Recovery"
%       IEEE Transactions on Signal Processing 51 nç7, (2003) p. 1771-1782
%
%   Copyright (c) 2006 Gabriel Peyre

options.null = 0;

if size(x0,1)==1 || size(x0,2)==1
    % 1D signal
    x0 = diff(x0); x0(end+1) = 0;
    x = perform_l1ball_projection(x0,tau);
    x = cumsum(x); x = x-mean(x)+mean(x0);
    err_tv = []; 
    err_l2 = [];
    return;
end

if isfield(options, 'mask')
    mask = options.mask;
else
    mask = [];
end

if isfield(options, 'x') && not(isempty(options.x))
    x = options.x;
else
    x = x0;
end
if isfield(options, 'niter')
    niter = options.niter;
else
    niter = 200;
end
err_tv = [];
err_l2 = [];
for i=1:niter
    progressbar(i,niter);
    % subgradient of the total variation
    t = perform_vf_normalization( grad(x, options) );
    t = -div( t, options );
    % gradient projection onto TV=tau
    tau1 = compute_total_variation(x, options);
    d = sum( t(:).^2 );
    if d>1e-9
        z = x - (tau1 - tau)*t / d;
    else
        z = x;
    end
    pi = dotp(x0-x,x-z);
    mu = dotp(x0-x,x0-x);
    nu = dotp(x-z,x-z);
    rho = mu*nu-pi^2;
    if rho==0 && pi>=0
        x = z;
    elseif rho>0 && pi*nu>=rho
        x = x0 + (1+pi/nu)*(z-x);
    elseif rho>0 && pi*nu<rho
        x = x + nu/rho*( pi*(x0-x)+mu*(z-x) );
    else
        error('PBM');
    end
    
    if not(isempty(mask))
        x = x.*mask;
    end
    
    % record errors
    err_tv(end+1) = tau1-tau;
    err_l2(end+1) = norm( x-x0 );
    if tau1<tau
        break;
    end
    
end


function vf = perform_vf_normalization(vf)

% perform_vf_normalization - renormalize a vf.
%
%   vf = perform_vf_normalization(vf);
%
%   Copyright (c) 2004 Gabriel Peyr?


d = sqrt( sum(vf.^2,3) );
d(d<1e-6) = 1;
vf = prod_vf_sf(vf,1./d);

function d = dotp(x,y)
d = sum(x(:).*y(:));

function v2 = prod_vf_sf(v1,s)

% prod_vf_sf - compute the product of a vector field by a scalar field.
%
%   v2 = prod_vf_sf(v1,s)
%
%   The result is the vector field defined by pointwise product.
%
%   Copyright (c) 2004 Gabriel Peyr?

v2 = v1 .* repmat(s, [1 1 2]);