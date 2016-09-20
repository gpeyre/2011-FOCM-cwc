% test for display of 1D eigenvectors and propagation

rep = ['results/display-1d/'];
if not(exist(rep))
    mkdir(rep);
end

lw = 2;
ratio = 4;

%% compute velocity field
n = 512;
x = (0:n-1)'/n;
freq = 2;
c = rescale( sin(2*pi*freq*x), 1,3);
[L,Delta] = compute_waveequation_matrix(c,'fft');

%% compute initial condition
options.vm = 0;
options.blurring_initial = 13; 
[f,df] = compute_initial_conditions(n,options);

%% display it
h = plot(c);
set(h, 'LineWidth', lw);
axis('tight'); set(gca, 'PlotBoxAspectRatio', [ratio 1 1]);
saveas(gcf, [rep 'speed-1d-' num2str(i) '.png'], 'png');

%% discretized time for resolution
tMax = cumsum(1./sqrt(c))/n; tMax = .5*max(tMax); % so that the whole domain is traveled once by the front
nt = 10; % number of timesteps
t = linspace(0,tMax, nt);
options.nplot = 6;
options.t = t;
options.tMax = tMax;
options.df = df;
options.tlist = t;

%% computing eigenvectors
[V,lambda] = svd(-diag(c)*Delta*diag(c));
lambda = diag(lambda);
[lambda,I] = sort(lambda,'ascend');
V = V(:,I);

%% display a few eigenvmodes
ne = 10;
elist = 2 + round( linspace(0,n/10,ne) );
a = max(max(abs(V(:,elist))));
for i=1:ne
    v = V(:,elist(i));
    clf;
    h = plot(v); 
    set(h, 'LineWidth', lw);
    axis([1 n -a a]); set(gca, 'PlotBoxAspectRatio', [ratio 1 1]);
    saveas(gcf, [rep 'eigenmode-1d-' num2str(i) '.png'], 'png');
end

%% solve for the evolution in the eigenvector domain
omega = sqrt(abs(lambda));
fw = V'*(f./c);
dfw = V'*(df./c);
Fw = repmat(fw,[1 nt]) .* cos( omega*t ) + repmat(dfw./omega,[1 nt]) .* sin( omega*t );
cT = repmat(c, [1 nt]);
F0 = cT .*(V*Fw);


%% display the solution
a = min(min(F0(:,2:end)));
b = max(max(F0(:,2:end)));
for i=1:nt
    v = F0(:,i);
    clf;
    h = plot(v); 
    set(h, 'LineWidth', lw);
    axis([1 n a b]); 
    set(gca, 'PlotBoxAspectRatio', [ratio 1 1]);
    saveas(gcf, [rep 'solution-1d-' num2str(i) '.png'], 'png');
end