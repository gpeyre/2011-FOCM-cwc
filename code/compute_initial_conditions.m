function [f,df] = compute_initial_conditions(n,options)

% compute_initial_conditions - compute initial conditions
%
%   [f,df] = compute_initial_conditions(n,options);
%
%   Copyrigth (c) 2007 Gabriel Peyre

ndim = 1;
if length(n)==2
    ndim = 2;
end

options.null = 0;
sigma = getoptions(options, 'blurring_initial', 7); % size of the bump, in pixels
vm = getoptions(options, 'vm', 0);

% position of the bump
c = getoptions(options, 'position_initial', 0.5);
c = round( (c-1/2)*n );

x = (-n/2+1:n/2)' - c(1);

% to convert into pixels
sigma = sigma/4;

if ndim==1
    if sigma==0
        f = zeros(n,1); f(round(n/2)) = 1;
    else
        f = exp( -x.^2/(2*sigma^2) );
    end
	df = zeros(n,1);
    if vm>0
        f = - (-1)^vm * diff(f, vm);
        f(end+1:n) = 0; 
    end
	f = f/max(f);
else
    y = (-n/2+1:n/2)' - c(2);
    if sigma==0
        f = zeros(n); 
        f(round(n/2),round(n/2)) = 1;
    else
        % initial value of the solution, a little gaussian bump
        [Y,X] = meshgrid(y,x);
        f = exp( -X.^2/(2*sigma^2)-Y.^2/(2*sigma^2) );
        if vm>0
            f = - f .* (X.^2+Y.^2-2*sigma^2);
            f = f-mean(f(:));
        end
    end
    f = f/max(f(:));
	% initial value of the derivative
	df = zeros(n);
end