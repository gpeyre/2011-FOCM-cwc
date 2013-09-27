% test for reverse time migration

path(path,'toolbox/');

if not(exist('n'))
    n = 1024;
    n = 2048/2;
end

name = 'twodiracs';

% width for display
lw = 2;
% font size
fs = 20;

rep = 'results/rtm/';
if not(exist(rep))
    mkdir(rep);
end


%% smooth medium

% contrast of the smooth medium
ampl = .05;
ctr = 1.4;
ampl = (ctr-1)/(ctr+1);
x = linspace(0,1,n+1)'; x(end) = [];
c0 = 1 + ampl*sin(2*pi*x);

%% medium + 2 bumps
dist = .08; % distance between spikes
b = [.5-dist/2 .5+dist/2];
h = [-1 1]*.5;  % height of the bumps
h = [-.5 .7];  % height of the bumps
sigma = .002;
c = c0;
for i=1:2
    g = exp( -(x-b(i)).^2/(2*sigma^2) );
    c = c + g*h(i);
end

%% for zooming for the display
Izoom = find( x>b(1)-.1 & x<b(2)+.1 );

clf;
h = plot(x, c); axis tight;
axis([0 1 0.3 1.6]);
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
saveas(gcf, [rep name '-rtm-speed-c.eps'], 'epsc');


% initial conditions
b0 = .03;
b0 = .01;

if not(exist('rtm_vm'))
    rtm_vm = 2;
end
options.vm = rtm_vm;
options.blurring_initial = 13; 
options.blurring_initial = 7; 
options.position_initial = b0;
[u0,du0] = compute_initial_conditions(n,options);

% u1 = -c0.*(u0([2:end,1])-u0([end 1:end-1]))*n/2;
u1 = -c0.*compute_spectral_derivative(u0);

derivative = 'fft';
[L,Delta] = compute_waveequation_matrix(c,derivative);


% discretized time for resolution
tMax = cumsum(1./c)/n; tMax = .7*max(tMax); % so that the whole domain is traveled once by the front
nt = 400; % number of timesteps
t = linspace(0,tMax, nt);

options.nplot = 6;
options.t = t;
options.tMax = tMax;
options.df = u1;
options.tlist = t;

% forward solution
disp('Computing fwd solution in c.');
F = compute_wave_solution(Delta,c, u0,u1, t);
plot_wave_solution(c, F, options);
saveas(gcf, [rep name '-rtm-fwd-c.eps'], 'epsc');

% isolate reflected wave
q0 = F(:,end);
a = .4;
q0(round(a*n):n) = 0;
q1 = compute_spectral_derivative(q0);
clf;
% subplot(2,1,1);
h = plot(x,q0); axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
% title('q0');
saveas(gcf, [rep name '-rtm-reflec-wave.eps'], 'epsc');
% subplot(2,1,2);
h = plot(x,q1); axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
title('q1');
saveas(gcf, [rep name '-rtm-reflec-wave-der.eps'], 'epsc');


% forward solution in smooth medium
disp('Computing fwd solution in c0.');
[F0,F0w,V,lambda] = compute_wave_solution(Delta,c0, u0,u1, t);
plot_wave_solution(c0, F0, options);
saveas(gcf, [rep name '-rtm-fwd-c0.eps'], 'epsc');

% backward solution in smooth medium
disp('Computing bwd solution in c0.');
[F1,F1w] = compute_wave_solution(Delta, c0, q0.*c,-q1.*c, t);
F1 = F1(:,end:-1:1);
plot_wave_solution(c0, F1, options);
saveas(gcf, [rep name '-rtm-bwd-c0.eps'], 'epsc');

% do the RTM
% 2/c^2(x)  *  int_0^T [dq/dt](x,t) [d^2u(/dt^2](x,t) dt
% backward*fwd''
disp('Perform RTM.');
R = perform_rtm_resolution(F0, F1, c0);

clf;
% Rp = -cumsum(R); Rp = Rp/max(Rp)*max(abs(c-c0));
h = plot(x, R); axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
saveas(gcf, [rep name '-rtm-true.eps'], 'epsc');


%% display zoom
clf;
h = plot(x(Izoom), c(Izoom)-c0(Izoom)); axis tight;
axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
saveas(gcf, [rep name '-rtm-residual-zoom.eps'], 'epsc');

Rmod = -cumsum(R)/n;
clf;
% Rp = -cumsum(R); Rp = Rp/max(Rp)*max(abs(c-c0));
h = plot(x(Izoom), Rmod(Izoom)); axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
saveas(gcf, [rep name '-rtm-zoom.eps'], 'epsc');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CWC RTM

nsub = 13;
nrandomization = 10;
options.nitermca = 500*4;
sublist = linspace(.8,.02,nsub);
options.thresh = 'hard';
options.thresh = 'soft';
options.solver = 'dr';

save_iter = 0;
err = zeros(nsub,nrandomization);
for k = 1:nsub
    % number of kept coefficients
    sub = sublist(k);
    m = round(sub*n);

    rand('seed', k*12345); % always same initia
    for irand=1:nrandomization
        if nrandomization>1
            progressbar(irand,nrandomization);
        end
        I = randperm(n); I = sort(I(1:m));
        % compressive basis
        Phi = V(:,I)';          % sub-set of eigenvectors
        lambdaI = lambda(I);    % sub-set of eigenvalues
      
        %% MCA computation
        % [F0,F0w,V,lambda] = compute_wave_solution(Delta,c0, u0,u1, t);
        options.F0 = F0;
        [F0T, err1] = perform_mca_cwc(Phi,lambdaI, u0,u1,            c0, t, options);
        
        % [F1,F1w] = compute_wave_solution(Delta,c0, q0.*c,-q1.*c, t);
        options.F0 = F1(:,end:-1:1);
        [F1T, err1] = perform_mca_cwc(Phi,lambdaI, q0.*c,-q1.*c,     c0, t, options);
        F1T = F1T(:,end:-1:1);
        
        RT = perform_rtm_resolution(F0T, F1T, c0);
        
        err(k,irand) = mean( (R(:)-RT(:)).^2 ) / mean( R(:).^2 );
    end
    
    if save_iter
        options.savefile = [repiter name '-rtm-sub' num2string_fixeddigit(round(100*sub),2)];
    end
    
    
%    RTp = -cumsum(RT); RTp = RTp/max(RTp)*max(abs(c-c0));
    clf;
    h = plot(RT); axis tight;
    set(h, 'LineWidth', lw);
    saveas(gcf, [rep name '-rtm-sub' num2string_fixeddigit(round(100*sub),2) '.eps'], 'epsc');
end

if 0
%% error averged accross realizations
errm = sqrt( mean(err, 2) ) / sqrt( mean(err(:).^2) );
clf;
h = plot(sublist,errm); axis tight;
set(h, 'LineWidth', lw);
saveas(gcf, [rep name '-rtm-error.eps'], 'epsc');
clf;
h = plot(sublist,log10(errm)); axis tight;
set(h, 'LineWidth', lw);
saveas(gcf, [rep name '-rtm-error-log.eps'], 'epsc');
end

