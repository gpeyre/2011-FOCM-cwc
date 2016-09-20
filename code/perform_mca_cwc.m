function [F, err] = perform_mca_cwc(Phi, lambdaI, f, df, c, t, options)

% perform_mca_cwc - MCA for compressive wave computations
%
%   [F, err] = perform_mca_cwc(Phi,lambdaI, f, df, c, t, options);
%
%   Phi is the matrix of a few eigenvectors
%   LambdaI is the corresponding eigenvalues
%   f and df are initial conditions at t=0
%   c is the speed of the medium
%   t is the time step
%
%   F is the solution at time t
%   err(i) it the error |F0-F| at iteration i of the algorithm
%
%   Copyright (c) 2008 Gabriel Peyre

t = t(:)';
nt = length(t);
n = length(c);
cT = repmat(c, [1 nt]);


% transform domain coefficients
fw = Phi*(f./c);    % fourier-like transform of the initial values
dfw = Phi*(df./c);  % fourier-like transform of the initial velocities
 
F0 = getoptions(options, 'F0', []);
thresh = getoptions(options, 'thresh', 'soft');
niter = getoptions(options, 'nitermca', 500);
solver = getoptions(options, 'solver', 'dr');

% solve for the evolution in the eigenvector domain, ie. compute <c^(-1)*f_t,phi_i>
omega = sqrt(abs(lambdaI));
Fw = repmat(fw,[1 nt]) .* cos( omega*t ) + repmat(dfw./omega,[1 nt]) .* sin( omega*t );

% for the remaining part, replace Phi by PhiW
W = getoptions(options, 'W', []);
PhiW = Phi;
if not(isempty(W))
    PhiW = Phi*W';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DOUGLAS RASHFORD CODE 

if strcmp(solver, 'dr')
    % solve for
    %  min |g|_1 subj. Phi*W'*g=Fw
    % and the solution is f=c*W'*g 
    options.pseudo_inv = PhiW';
    options.x0 = F0;
    options.verb = 0;
    options.niter = niter;
    [F,lun,Err] = perform_douglas_rachford(PhiW, Fw, options);
    if not(isempty(W))
        F = W'*F;
    end
    F = F .* cT;
    % compute error
    err = mean( (F0(:)-F(:)).^2 ) /  mean(F0(:).^2); 
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MCA CODE 

if 0
    % back to the spacial domain (L^2 reconstruction)
    F = pinv(Phi)*Fw;
    % solution with L2 unregularized inverse
    plot_wave_solution(c,F,options);
    saveas(gcf, [rep name  '-rough' num2str(round(100*roughness)) '-nsub' num2str(k) '-pinv.png'], 'png');
end

use_mad = getoptions(options, 'use_mad', 0);

% compute the set of thresholds
if 0
    % use simple max
    tmax = max( max( abs(cT.*(Phi'*Fw)) ) ); % max of thresholds
    tlist = linspace(tmax, 0, niter);
else
    % use more clever median
    sigma_mca = 3;
    x = cT.*(Phi'*Fw);
    tmax = mad( x(:), 1 )/0.6745 * sigma_mca;
%    tlist = tmax' * linspace(1, 0, niter);
end

tlist = linspace(tmax, 0, niter);

% solve for min |c.*g|  subj.t. Phi*g=Fw,   and the solution is f_t=c*g
mu = 1; % step size 
G1 = zeros(n,nt);   % solution at various time t
err = []; % see in oracle how it decays
for i=1:niter
    % project
    G1 = G1 + Phi'*( Fw-Phi*G1 );
    % set up threshold
    tval = tlist(i);
    
    if use_mad
        x = G1.*cT;
        tval = mad( x(:), 1 )/0.6745 * sigma_mca;
    end
    
    % threshold
    G1 = perform_thresholding(G1.*cT,tval,thresh)./cT;
    F1 = G1.*cT;
    % compute error
    if not(isempty(F0))
        err(end+1) = mean( (F0(:)-F1(:)).^2 ) /  mean(F0(:).^2); % / norm(F0,'fro');
        if i==1 || err(i)<min(err(1:i-1))
            F = F1; % record minium error solution
        end
    else
        F = F1;
    end
end