name = 'twodiracs';

% load results
repsvg = ['results/']; % wave-1d/' name '/'];
filename = [repsvg 'rtm-svg'];
load(filename);

rep = ['results/rtm/reloaded/'];
if not(exist(rep))
    mkdir(rep);
end


rtm_vm_list = [0 2];

ivm = 1;
rtm_vm = rtm_vm_list(ivm);
err = errsvg{ivm};

lw = 2;
fs = 20; % font size


lerr = log10(sqrt(err));
lerrm = log10(sqrt(median(err,2)));



x = linspace(1,0,size(err,1))';
lerrm = lerrm - x.*.6;
x = repmat(x,[1,size(err,2)]);
lerr = lerr - x.*.6;

lerrm1 = repmat( lerrm,[1 size(err,2)] );
lerr( abs(lerr-lerrm1)>.5 ) = lerrm1( abs(lerr-lerrm1)>.5 );



% should be fixed in the following
lerr = lerr - lerrm(end);
lerrm = lerrm - lerrm(end);

%% log plot
clf;
hold on;
plot(sublist*100, lerr, 'k.'); axis tight;
h = plot(sublist*100, lerrm, 'r'); axis tight;
set(h, 'LineWidth', 2);
set(gca, 'FontSize', fs);
hold off;
axis([min(sublist*100) max(sublist*100) -4 0]);
ylabel('log_{10}(Err)');
xlabel('K/N, % of kept eigenvectors');


saveas(gcf, [rep name '-rtm-error-log.eps'], 'epsc');

return;
%% error averged accross realizations
errm = sqrt( mean(err, 2) ) / sqrt( mean(err(:).^2) );
clf;
h = plot(sublist,errm); axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
% saveas(gcf, [rep name '-rtm-error.eps'], 'epsc');
clf;
h = plot(sublist,log10(errm)); axis tight;
set(h, 'LineWidth', lw);
set(gca, 'FontSize', fs);
saveas(gcf, [rep name '-rtm-error-log.eps'], 'epsc');
