function M = compute_transform_matrix(w, callb, options)

% compute_haar_matrix - compute the matrix for an orthogonal transform given using a matrix.
%
%	M = compute_haar_matrix(w, callb, options);
%
% if options.ndims = 2 : 2D transform
%	w is the size of the square (so that M is a transform matrix
%	of size w^2 x w^2).
% if options.ndims = 1 (default) : 1D transform
%
%   The forward transform is computed via y=M*x
%   so each row of M is a basis atom.
%
%   Copyright (c) 2008 Gabriel Peyrée

options.null = 0;

ndims = getoptions(options, 'ndims', length(w));
verb = getoptions(options, 'verb', 1);

if ndims==1
    M = zeros(w,w);
    for i=1:w
        y = dirac(w, i);
        x = feval(callb, y, options);
        M(:,i) = x(:);
    end
    M = M';
    return;
end

if length(w)==1
    w= [w w];
end

n = prod(w);
M = zeros(n,n);
k = 0;
for i=1:w(1)
    for j=1:w(2)
        k = k+1;
        if verb
            progressbar(k,n);
        end
        y = dirac(w, [i,j]);
        x = feval(callb, y, options);
        M(k,:) = x(:)';
    end
end

function y = dirac(s,n)

% dirac - dirac function.
%
%   y = dirac(s,n);
%
%   s is the size of the matrix, and 
%   n is the location of the dirac (default (1,...,1)).
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    n = ones(size(s));
end

if length(s)==1
    y = zeros(s,1);
else
    y = zeros(s);
end
y = array_set_val( y, n, 1 );

function y = array_set_val( x, ind, v )

% array_set_val - set a value in a multidimensional array using a vector for the index
%
%   y = array_set_val( x, ind, v );
%
%   Copyright (c) 2004 Gabriel Peyré

ind = num2cell(ind);
x( ind{:} ) = v;
y = x;

