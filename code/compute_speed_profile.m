function c = compute_speed_profile(name,n,options)

% compute_speed_profile - load a speed domain
%
%   c = compute_speed_profile(name,n,options);
%
%   n=[n1 n2] for 2D domains (otherwise 1D).
%
%   A blurring of width options.blurring (in pixels) is performed as a post processing.
%
%   Two parameters controls the complexity of the speed function:
%   options.roughness should be in [1,10], controls complexity of the medium
%   options.contrast should be >1 and controls max/min value of the medium
%
%   Copyright (c) 2007 Gabriel Peyre

ndim = 1;
if length(n)==2
    ndim = 2;
end

options.null = 0;
sigma       = getoptions(options, 'sigma', .06);
blurring    = getoptions(options, 'blurring', 15);
roughness    = getoptions(options, 'roughness', 1);
contrast    = getoptions(options, 'contrast', 2);

if ndim==1
    %%%%% 1D %%%%%
    x1 = linspace(-1,1,round(n/3))';
    x2 = linspace(-1,1,round(n/3))';
    x3 = linspace(-1,1,n-2*round(n/3))';
    switch lower(name)
        case 'constant'
            c = ones(n,1);
        case {'sin' 'sine'}
            x = (0:n-1)/n;
            f = 2*round(roughness); % between 1/10 oscillations
            c = sin(2*pi*f*x);
        case 'piece-constant'
            c = [x1*0+1; x2*0+2; x3*0+3];
%            c = [x1*0+1; x2*0+1; x3*0+3];
        case 'piece-polynomial'
            c = [x1.^2; 1/2 + x2/2; 1/2 + (x3.^3)/2];       
        case 'rand'
            c = randn(n,1)*sigma+1;
        case 'bv'
            tv_divide = max(1, 50*(1-roughness) );
            c = randn(n,1);
            tv = compute_total_variation(c, options);
            tau = tv/tv_divide;
            options.bound = 'per';
            [c,err_tv,err_l2] = perform_tv_projection(c,tau,options);
            tv1 = compute_total_variation(c, options);
        case {'regsteps' 'steps'}
            % number of steps, should be even
            nsteps = 2*roughness; % between 2 and 20 steps
            x = (0:n-1)/n + 1/(2*nsteps);
            c = sin(nsteps*pi*x)>=0;
        case 'steps-varying'
            x = (0:n-1)'/n;
            % number of steps
            m = 2*ceil(roughness);
            % signs
            s = (-1).^((1:m)');
            % locations: uniform
            t = round( linspace(1,n,m+1) + n/(2*m) ); t(end) = [];
            % location: varying
            t = linspace(0,1,m+1).^2 + .5*linspace(0,1,m+1); 
            t = rescale(t(:),1,n);  t = round(t+t(2)/2); t(end) = [];            
            % speed without modulation
            c = zeros(n,1); c(t) =  s;
            c = rescale(cumsum(c),-1,1);
            % modulate
            cm = rescale( sin(pi*x + pi/8).^2, .5, 1);            
            c = c.*cm;
        case 'bv-steps'
            % number of steps
            m = 1 + floor(roughness);
            % steps values
            s = 1 + .2*randn(m,1); % rand(m,1)*2-1;
            s = (-1).^((1:m)') .* s;
            % change of sign to ensure in [-1,1]
            if 0
            step_bounds    = getoptions(options, 'step_bounds', 1);
            v = 0;
            for i=1:m
                if abs(v+s(i))>step_bounds
                    s(i) = -s(i);
                end
                v = v+s(i);
            end
            end
            % set random position
            dx = rescale(rand(m+2,1),.3,1);
            q = cumsum(dx); q = rescale(q,1,n); q(1)=[]; q(end)=[];           
%            q = randperm(n); q = sort( q(1:m) );
            c = zeros(n,1); c(round(q)) = s;
            c = cumsum(c);
        otherwise
            error('Unknown speed name');
    end
    c = c(:);
else
    %%%% 2D %%%%%
    switch lower(name)
        case 'constant'
            c = ones(n);
        case 'sin'
            x = (0:n(1)-1)/n(1);
            y = (0:n(2)-1)/n(2);
            [Y,X] = meshgrid(y,x);
            f = round(roughness);
            c = sin(2*pi*f*X) .* sin(2*pi*f*Y);
        case 'layered'
            c = randn(n,1)*sigma+1;
            c = repmat(c, 1, n);
        case 'layered2'
            c = max( 0.05, randn(n,1)*4*sigma+1 );
            c = repmat(c, 1, n);
        case 'random'
            c = rand(n)*sigma + 1;
        case 'bv'
            tv_divide = max(1, 50*(1-roughness) );
            c = randn(n);
            tv = compute_total_variation(c, options);
            tau = tv/tv_divide;
            options.niter = min(150*tv_divide,3000);
            options.bound = 'per';
            [c,err_tv,err_l2] = perform_tv_projection(c,tau,options);
            tv1 = compute_total_variation(c, options);
            disp( ['... tv error: ' num2str( abs(tv1-tau)/tau*100, 2 ) '%' ] );
            x = linspace(1,contrast,prod(n));
            c = perform_histogram_equalization(c,x);
            % c = clamp(1+.9*(c-mean(c(:)))/(2*std(c(:))),1-eta,1+eta);
        case 'disk'
            x = linspace(0,1,n(1));
            [Y,X] = meshgrid(x,x);
            c = [.4,.6]; r = .3;
            c = (X-c(1)).^2 + (Y-c(2)).^2 < r^2;
    end
end

% post blurring
if blurring>0
    options.bound = 'per';
    c = perform_blurring(c,blurring, options);
end

c = rescale(c, 1, contrast);