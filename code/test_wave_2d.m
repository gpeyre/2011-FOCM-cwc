% test for wave equation in 2D using compressive sensing in the eigenbasis

path(path,'toolbox/');

% size of the problem
if not(exist('n'))
    n = 32*2;
end

% speed function
if not(exist('name'))
name = 'constant';
name = 'layered2';
name = 'layered';
name = 'bv';
name = 'disk';
name = 'sin';
end

rep = ['results/wave-2d/' name '/'];
if not(exist(rep))
    mkdir(rep);
end
repiter = ['results/wave-2d/' name '/iter/'];
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
    speed_blur = 2;
end
if not(exist('blurring_initial'))
    blurring_initial = 4;
end
if not(exist('thresh'))
    thresh = 'soft';
end
if not(exist('nrandomization'))
    nrandomization = 1;
end
if not(exist('nsub')) % number of sub-sampling for #eigenvectors
    nsub = 10;
end
if not(exist('nitermca'))
    nitermca = 200;
end
if not(exist('save_iter'))  % save progressions
    save_iter = 1;
end
if not(exist('add_str'))
    add_str = '';
end
if not(exist('save_eps'))
    save_eps = 0;
end
if not(exist('additional_plot'))
    additional_plot = 1;
end
if not(exist('domain'))
    % sparsity transform
    domain = 'spac';
    domain = 'wav';
end
if not(exist('eig_sampling'))
    eig_sampling = 'uniform';
    if strcmp(domain, 'wav')
        eig_sampling = 'biased';
    end
end

options.roughness = roughness;
options.contrast = contrast;
options.blurring = speed_blur;  % in pixels
disp('--> Compute speed.');
c = compute_speed_profile(name,[n n],options);

% colormap for the animation
options.colormap = jet(256);

% initial value of the solution, a little gaussian bump
options.blurring_initial = blurring_initial; % size of the bump
options.vm = 2; %  use laplacian of Gaussian
[f,df] = compute_initial_conditions([n n],options);

% to name files
str = ['-rough' num2str(roughness) '-contrast' num2str( round(contrast*10) ) add_str];


% derivative operator
derivative = 'fft';
disp('--> Computing wave matrix.');
[L,Delta] = compute_waveequation_matrix(c,derivative);

% discretized time for resolution
tMax = cumsum(1./sqrt(c(:,end/2)))/n; 
tMax = cumsum(1./(c(:,end/2)))/n; 
% tMax = .5*max(tMax); % so that the whole domain is traveled once by the front
tMax = .5*max(tMax); % so that the whole domain is traveled once by the front
if not(exist('nt'))
    nt = 40; % number of timesteps
end
t = linspace(0,tMax, nt);
options.nplot = 6;
options.t = t;

options.tMax = tMax;
options.df = df;
options.tlist = t;

% eigenbasis, L=-diag(c)*V*diag(lambda)*V'*diag(c.^(-1))
disp('--> Computing eigenvectors.');
[V,lambda] = svd(-diag(c(:))*Delta*diag(c(:)));
lambda = diag(lambda);

% put low frequencies first
[lambda,I] = sort(lambda, 'descend');
V = V(:,I);
    
% solve for the evolution in the eigenvector domain
% (this is an oracle solution)
disp('--> Computing inverse matrix');
omega = sqrt(abs(lambda));
fw = V'*(f(:)./c(:));
dfw = V'*(df(:)./c(:));
Fw = repmat(fw,[1 nt]) .* cos( omega*t ) + repmat(dfw./omega,[1 nt]) .* sin( omega*t );
cT = repmat(c(:), [1 nt]);
F0 = cT .*(V*Fw);
F0 = reshape(F0,[n n nt]);

% save movie
compute_movie_file( normalize_movie(F0) , [rep name '-true.avi'], options);
% save images
plot_wave_solution(c, normalize_movie(F0) ,options);
% saveas(gcf, [rep name '-true.png'], 'png');

% progression of sub-sampling
if not(exist('sublist')) || length(sublist)~=nsub
    sublist = linspace(.8,.02,nsub);
end
nsub = length(sublist);

% thresholding operator used for MCA iterations
options.thresh = thresh;

%% compute explicitely the matrix of the transformation
if strcmp(domain,'spac')
    W = eye(n^2,n^2);
else
    disp('--> Computing transform matrix');
    % options for the wavelet transform
    options.wavelet_type = 'daubechies';
    options.wavelet_vm = 4;
    Jmin = 1;
    callb = @(x,options)perform_wavelet_transform(x,1,-1,options);
    W = compute_transform_matrix([n n], callb, options);
end
options.W = W;

err = zeros(nsub,nrandomization);
err_lin = zeros(nsub,nrandomization);
for k = 1:nsub
    % number of kept coefficients
    sub = sublist(k);
    m = round(sub*n^2);
    rand('seed', k*12345); % always same initialization
    for irand=1:nrandomization
        if nrandomization>1
            progressbar(irand,nrandomization);
        end        
        if strcmp(eig_sampling, 'uniform')       
            I = randperm(n^2); I = sort(I(1:m));
        else
            J = 6;  % number of constant cluster
            I = [];
            N = n^2; 
            L = N/2^(J-1); D = 0;
            for j=1:J
                sel = D+1:D+L;
                mj = min(floor(m/J),length(sel));
                if j==J
                    mj = m - length(I);
                    sel = D+1:N;
                end
                I1 = randperm(length(sel)); I1 = sort( sel(I1(1:mj)) );
                I = [I; I1(:)];
                % advance
                D = D+L;
                if j>1
                    L = L*2;
                end
            end
            I = sort(I);
        end
        
        % compressive basis
        Phi = V(:,I)';          % sub-set of eigenvectors
        lambdaI = lambda(I);    % sub-set of eigenvalues
        
        %% MCA/DR computation over the transformed domain
        options.F0 = reshape( F0, [n*n nt]);
        options.nitermca = nitermca; 
        [F1, err1] = perform_mca_cwc(Phi, lambdaI, f(:), df(:), c(:), t, options);
        F1 = reshape(F1, [n n nt]);
        % save for later display
        Fsvg{k} = F1;        
        % the error is the squared error
        err(k,irand) = min(err1);        
    end
    
    % save solution
    str = ['-rough' num2str(round(100*roughness)) '-sub' num2string_fixeddigit(round(100*sub),2)];  
    str = ['-sub' num2string_fixeddigit(round(100*sub),2) '-' domain];
    plot_wave_solution(c, normalize_movie(Fsvg{k}) ,options);
    % saveas(gcf, [rep name str '.png'], 'png');
    compute_movie_file( normalize_movie(Fsvg{k}) , [rep name str '.avi'], options);
end


if exist('additional_plot') && additional_plot    
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
    
    saveas(gcf, [rep name '-error-log-' domain '.png'], 'png');
    if save_eps
        saveas(gcf, [rep name str '-error-log-' domain '.eps'], 'epsc');
    end
end

