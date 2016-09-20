% test for 1D wave equation using compressive sensing in the eigenbasis
%
% d^2(f_t)/dt^2 = L*f_t    f_0=f
%   L=c^2*Delta, c diagonal matrix of speed
% Factor c*Delta*c = -V*S*V', V ortho-eigenvector, S diagonal values
%   L=c^2*Delta = -c*V*S*V'*c^{-1}
% => d^2(g_t)/dt^2=-S*g_t  with  g_t=V'*c^{-1}*f_t
%   V=[v_i]_{i},   <c^{-1}*f_t,v_i> = <c^{-1}*f,v_i>*cos( sqrt(s(i))*t )
%
% from a set of M computed values 
%      y=[<c^{-1}*f_t,v_{a(i)}>]_{i=1}^M=Phi*f_t
%    where Phi=[v_{a(i)}]_i'  (v_{a(i)} are rows)
% on recover f_t=c*g_t by the following optimization
%   min |c*g_t|_1   subj.to   Phi*g_t=y
%
%

path(path,'toolbox/');

% size of the problem
if not(exist('n'))
    n = 256*2;
    n = 512*2;
    n = 2048;
end

% speed function
if not(exist('name'))
    name = 'constant';
    name = 'rand';
    name = 'piece-polynomial';
    name = 'bv';
    name = 'piece-constant';
    name = 'regsteps';
    name = 'bv-steps';
    name = 'sin';
end

rep = ['results/wave-1d/' name '/'];
if not(exist(rep))
    mkdir(rep);
end
repiter = ['results/wave-1d/' name '/iter/'];
if not(exist(repiter))
    mkdir(repiter);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% options
if not(exist('roughness'))
    roughness = 5;
end
if not(exist('contrast')==1)
    contrast = 3;
end
if not(exist('speed_blur'))
    speed_blur = 5;
end
if not(exist('blurring_initial'))
    blurring_initial = 11;
end
if not(exist('thresh'))
    thresh = 'soft';
end
if not(exist('nrandomization'))
    nrandomization = 10;
end
if not(exist('nsub')) % number of sub-sampling for #eigenvectors
    nsub = 15;
end
if not(exist('nitermca'))
    nitermca = 300;
end
if not(exist('save_iter'))  % save progressions
    save_iter = 0;
end
if not(exist('add_str'))
    add_str = '';
end
if not(exist('save_eps'))
    save_eps = 0;
end

options.roughness = roughness;
options.contrast = contrast;
options.blurring = speed_blur;  % in pixels
disp('--> Compute speed.');
c = compute_speed_profile(name,n,options);


% to name files
str = ['-rough' num2str(roughness) '-contrast' num2str( round(contrast*10) ) add_str];

% derivative operator
derivative = 'fd';
derivative = 'fft';
disp('--> Computing wave matrix.');
[L,Delta] = compute_waveequation_matrix(c,derivative);

% initial value of the solution, a little gaussian bump
options.blurring_initial = blurring_initial; % size of the bump
options.vm = 0;
if strcmp(name, 'sin')
    options.vm = 2;
    options.blurring_initial = 13; 
end
[f,df] = compute_initial_conditions(n,options);

% discretized time for resolution
tMax = cumsum(1./sqrt(c))/n; tMax = .5*max(tMax); % so that the whole domain is traveled once by the front
nt = 20; % number of timesteps
t = linspace(0,tMax, nt);
options.nplot = 6;
options.t = t;

options.tMax = tMax;
options.df = df;
options.tlist = t;

if 0
    disp('--> Solving with RK4.');
    [Frk,trk] = compute_waveequation_rk(L, f, options);
end

% eigenbasis, L=-diag(c)*V*diag(lambda)*V'*diag(c.^(-1))
disp('--> Computing eigenvectors.');
if 0
    [V,lambda] = eig(-diag(c)*Delta*diag(c));
    lambda = real(diag(lambda)); % eigenvalues
    % lambda = diag(lambda);
    % impose real and orthogonality
    [V,R] = qr(real(V));
else
    [V,lambda] = svd(-diag(c)*Delta*diag(c));
    lambda = diag(lambda);
end

% solve for the evolution in the eigenvector domain
% (this is an oracle solution)
disp('--> Computing inverse matrix');
omega = sqrt(abs(lambda));
fw = V'*(f./c);
dfw = V'*(df./c);
Fw = repmat(fw,[1 nt]) .* cos( omega*t ) + repmat(dfw./omega,[1 nt]) .* sin( omega*t );
cT = repmat(c, [1 nt]);
F0 = cT .*(V*Fw);

options.savefile = '';
if save_iter
    options.savefile = [repiter name str '-true'];
    plot_wave_solution(c, F0, options);
end
if isempty(options.savefile) && save_iter
    saveas(gcf, [repiter name str '-true.png'], 'png');
end
    
% progression of sub-sampling
if not(exist('sublist')) || length(sublist)~=nsub
    sublist = linspace(.8,.02,nsub);
end
%    sublist = [.1 .2 .3 .4 .5];
% nsub = length(sublist);

% thresholding operator used for MCA iterations
options.thresh = thresh;

err = zeros(nsub,nrandomization);
err_lin = zeros(nsub,nrandomization);
for k = 1:nsub
    % number of kept coefficients
    sub = sublist(k);
    m = round(sub*n);

    rand('seed', k*12345); % always same initialization
    for irand=1:nrandomization
        if nrandomization>1
            progressbar(irand,nrandomization);
        end
        I = randperm(n); I = sort(I(1:m));
        % compressive basis
        Phi = V(:,I)';          % sub-set of eigenvectors
        lambdaI = lambda(I);    % sub-set of eigenvalues
        
        %% MCA computation
        options.F0 = F0;
        options.nitermca = nitermca;
        [F1, err1] = perform_mca_cwc(Phi,lambdaI, f, df, c, t, options);

        % the error is the squared error
        Fsvg{k} = F1;
        err(k,irand) = min(err1);        
    end
    
    if save_iter
        options.savefile = [repiter name str '-sub' num2string_fixeddigit(round(100*sub),2)];
        plot_wave_solution(c,Fsvg{k},options);
    end    
    if isempty(options.savefile) && save_iter
        saveas(gcf, [repiter name str '-sub' num2string_fixeddigit(round(100*sub),2) '.png'], 'png');
    end
end

%% regular plot
if 0
    clf;
    hold on;
    plot(sublist*100, err, 'k.'); axis tight;
    h = plot(sublist*100, mean(err,2), 'r'); axis tight;
    set(h, 'LineWidth', 2);
    hold off;
    axis([min(sublist*100) max(sublist*100) 0 1]);
    ylabel('Rel. L^2 error');
    xlabel('% of kept eigenvectors');
    saveas(gcf, [rep name str '-error.png'], 'png');
    if save_eps
        saveas(gcf, [rep name str '-error.eps'], 'epsc');
    end
end

if not(exist('additional_plot'))
    additional_plot = 0;
end

if additional_plot
    
    lerr = log10(sqrt(err));  
    lerrm = log10(sqrt(median(err,2)));
    lerr = lerr - lerrm(end);
    lerrm = lerrm - lerrm(end);

    %% log plot
    clf;
    hold on;
    plot(sublist*100, lerr, 'k.'); axis tight;
    h = plot(sublist*100, lerrm, 'r'); axis tight;
    set(h, 'LineWidth', 2);
    set(gca, 'FontSize', 20);
    hold off;
    axis([min(sublist*100) max(sublist*100) -5 0]);
    ylabel('log_{10}(Err)');
    xlabel('K/N, % of kept eigenvectors');
    
    saveas(gcf, [rep name str '-error-log.png'], 'png');
    if save_eps
        saveas(gcf, [rep name str '-error-log.eps'], 'epsc');
    end
end