function x = compute_rand_sparse(m,s, issigned)

% x = compute_rand_sparse(m,s, issigned);

if nargin<3
    issigned = 0;
end

% random +/- 1 signal
x = zeros(m,1);
q = randperm(m);
if issigned
    x(q(1:s)) = sign(randn(s,1));
else
    x(q(1:s)) = 1;
end