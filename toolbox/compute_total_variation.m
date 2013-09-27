function TV = compute_total_variation(y, options)

% compute_total_variation - compute the total variation of an image
%
%   TV = compute_total_variation(y);
%
%   Copyright (c) 2007 Gabriel Peyre

nbdims = 2;
if size(y,1)==1 || size(y,2)==1
    nbdims = 1;
end

if nbdims==1
    TV = sum( abs(diff(y)) );
    return;
end

% options.bound = 'per';
options.null = 0;
g = grad(y, options);
TV = sqrt( sum(g.^2,3) );
TV = sum(TV(:));