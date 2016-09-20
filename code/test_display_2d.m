% test for display of 2D eigenvectors and propagation

rep = ['results/display-2d/'];
if not(exist(rep))
    mkdir(rep);
end

lw = 2;
ratio = 4;

%% compute velocity field
n = 40;
x = (0:n-1)'/n;
[Y,X] = meshgrid(x,x);
freq = 1;
c = rescale( sin(2*pi*freq*X) .* sin(2*pi*freq*Y) , 1,3);
[L,Delta] = compute_waveequation_matrix(c,'fft');

%% compute initial condition
options.vm = 0;
options.blurring_initial = 4; 
[f,df] = compute_initial_conditions([n n],options);

%% display it
imageplot(c);
save_image( rescale(c), [rep 'speed-2d']);

%% discretized time for resolution
tMax = cumsum(1./(c(:,end/2)))/n; 
tMax = .5*max(tMax); % so that the whole domain is traveled once by the front
nt = 10; % number of timesteps
t = linspace(0,tMax, nt);
options.nplot = 6;
options.t = t;
options.tMax = tMax;
options.df = df;
options.tlist = t;

%% computing eigenvectors
[V,lambda] = svd(-diag(c(:))*Delta*diag(c(:)));
lambda = diag(lambda);
[lambda,I] = sort(lambda,'ascend');
V = V(:,I);

%% display a few eigenvmodes
ne = 10;
elist = 3 + round( linspace(0,n^2/10,ne) );
options.base_str = [rep 'eigenmode-2d-'];
for i=1:ne
    v = reshape( V(:,elist(i)), n,n);
    save_image( rescale(v), num2str(i), options);
end

%% solve for the evolution in the eigenvector domain
omega = sqrt(abs(lambda));
fw = V'*(f(:)./c(:));
dfw = V'*(df(:)./c(:));
Fw = repmat(fw,[1 nt]) .* cos( omega*t ) + repmat(dfw./omega,[1 nt]) .* sin( omega*t );
cT = repmat(c(:), [1 nt]);
F0 = cT .*(V*Fw);
F0 = reshape(F0,[n n nt]);

%% display the solution
options.base_str = [rep 'solution-2d-'];
a = min(F0(:)); b = max(F0(:));
for i=1:nt
    v = F0(:,:,i);
    save_image( rescale(v), num2str(i), options);
end