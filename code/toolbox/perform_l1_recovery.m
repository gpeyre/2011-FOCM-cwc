function x = perform_l1_recovery(D,y,options)

options.null = 0;

if isfield(options, 'method')
    method = options.method;
else
    method = 'bp';
end
if isfield(options, 's')
    s = options.s;
else
    s = 20;
end

m = size(D,2);

% OMP options
maxItersOMP = s;
solFreq = 0; verbose = 0; lambdaStop = 0;
% BP options
maxIters = 20; lambda = 0; OptTol = 1e-5;
% MCA options
mcaIters = 30;
if isfield(options, 'lambda')
    lambda = options.lambda;
end
if isfield(options, 'mcaIters')
    mcaIters = options.mcaIters;
end

switch lower(method)
    case 'bp'
        for i=1:size(y,2)
            x(:,i) = SolveBP(D, y(:,i), m, maxIters, lambda, OptTol);
        end
    case 'omp'
        % x = SolveOMP(D, y, m, maxItersOMP, lambdaStop, solFreq, verbose, OptTol);
        options.nbr_max_atoms = s;
        x = perform_omp(D,y,options);
    case 'mca'
        options.niter = mcaIters;
        [x,T] = perform_sparse_inversion(D,y,options);
    otherwise
        error('Unknown method');
end