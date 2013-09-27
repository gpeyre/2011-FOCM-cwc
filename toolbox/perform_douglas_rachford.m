function [x1,lun,err] = perform_douglas_rachford(A,y,options)

% perform_douglas_rachford - 
%
%   [x1,lun,err] = perform_douglas_rachford(A,y,options);
%
%   Solve noiseless Basis Pursuit
%       x1 = argmin_x \sum |x[k]|   subj.to A*x=y
%
%   A can be a matrix or an operator A(x,dir,options)
%   with dir=1 for A*x and dir=-2 for A^{+}*x (pseudo inverse).
%
%   Special thanks to Jalal Fadili for useful helps on this algorithm.
%
%   Copyright (c) 2008 Gabirel Peyre

niter = getoptions(options, 'niter', 100);
mu = getoptions(options, 'mu', 3);

ntest = size(y,2);

if isnumeric(A)
    n = size(A,2);
    pA = getoptions(options, 'pseudo_inv', []);
    if isempty(pA)
        pA = pinv(A);
    end
else
    % retrieve size
    x = feval( A, y*0, -1, options );
    n = size(x,1);
end

x = getoptions(options, 'x', zeros(n,ntest) );
x0 = getoptions(options, 'x0', zeros(n,ntest) );

drawiter = getoptions(options, 'drawiter',0);
verb = getoptions(options, 'verb', 1);

lun = [];
err = [];
if isnumeric(A)
    x1 = x + pA*(y-A*x);
else
    x1 = y - feval(A,x, +1, options);
    x1 = x + feval( A, x1, -1, options);
end
if drawiter
    clf;
end

nrefresh = 100;

for i=1:niter
    if verb
        progressbar(i,niter);
    end
    % compute the intermediary point
    %   x <- x-x1 + thresh( 2*x1-x )
    u = real(2*x1-x); % the vector to threshold
    s = abs(u) - mu; s = (s + abs(s))/2;
    u = sign(u).*s;    
    x = x-x1 + u;
    % perform projection step on the constraints
    if isnumeric(A)
        x1 = x + pA*(y-A*x);
    else
        x1 = y - feval(A,x,+1, options);
        x1 = x + feval( A, x1, -2, options);
    end
    % x1 = real(x1); % force real signal
    % check for error decay
    if nargout>=3
        err(end+1) = norm(x0-x1, 'fro')^2;
    end
    if nargout>=2
        lun(end+1) = mean(sum(abs(x1)));
    end
    if drawiter && mod(i,ceil(niter/nrefresh))==1  % && ntest==1
        plot( real(x1(:,1:min(end,3))) ); axis([1,n,-1,1]); drawnow;
    end
end
x1 = real(x1);
