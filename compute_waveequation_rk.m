function [F,t] = compute_waveequation_rk(L, f, options)

% compute_waveequation_rk - solve the wave equation with RK4 time stepping
%
%   [F,t] = compute_waveequation_rk(L,f,options);
%
%   Solve
%       d^2F/dt^2 = L*F    for t in [0,options.tMax]
%       with initial condition F(t=0)=f and dF/dt(t=0)=options.df
%   with RK4/5 solver.
%
%   F(:,k) is the solution at time t(k).
%
%   Copyright (c) Gabriel Peyre


tMax    = getoptions(options, 'tMax', .4);
df      = getoptions(options, 'df', f*0);
tlist   = getoptions(options, 'tlist', []);

% global L;
odefun = @(t,x) [x(end/2+1:end); L*x(1:end/2)];
Tspan = [0,tMax];
F0 = [df;f];
[t,F] = ode45(odefun,Tspan,F0);
F = F(:,end/2+1:end)';

if not(isempty(tlist))
    F1 = F;
    n  = size(F,1);
    nt = length(tlist);
    F = zeros(n,nt);
    for i=1:n
        F(i,:) = interp1(t,F1(i,:),tlist);
    end
end