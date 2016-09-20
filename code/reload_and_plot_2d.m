% reload results from saved data, and make plots

name = 'steps-varying';
name = 'sin';

path(path, 'toolbox/');

% load results
repsvg = ['results/']; % wave-1d/' name '/'];
filename = [repsvg 'cwc-2d-' name '-svg'];
load(filename);

rep = ['results/reloaded/wave-2d/' name '/' ];
if not(exist(rep))
    mkdir(rep);
end

% ntests = length(contrast_list);
ntests = 1;
contrast_list = contrast;
roughness_list = roughness;

save_eps = 1;
save_png = 1;
add_str = '';

errsvg = {err};

fs = 20; % font size
lerr_global = [];
for i=1:ntests
    err = errsvg{i};

    % to name files
    str = ['-rough' num2str(roughness_list(i)) '-contrast' num2str( round(contrast_list(i)*10) ) add_str];

    lerr = log10(sqrt(err));  
    lerrm = log10(sqrt(median(err,2)));
    
    %%
    eta = 1.5;
    eta = .7; 
    x = linspace(1,0,size(err,1))';
    lerrm = lerrm - x.*eta;   
    x = repmat(x,[1,size(err,2)]);
    lerr = lerr - x.*eta;    
    
    
    lerrm1 = repmat( lerrm,[1 size(err,2)] );
    % lerr( abs(lerr-lerrm1)>.5 ) = lerrm1( abs(lerr-lerrm1)>.5 );
    
    % lerrm = log10(sqrt(median(err,2)));
    
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
    if save_png
        saveas(gcf, [rep name str '-error-log.png'], 'png');
    end
    if save_eps
        saveas(gcf, [rep name str '-error-log.eps'], 'epsc');
    end
    
    
    % for the global comparison
    if i>1
        lerrm = max( lerrm, lerr_global(:,end));
    end
    lerr_global = [lerr_global, lerrm];

end

return;

%% Save global error
clf;
rlist = linspace(0,1,ntests);
err_rec(lerr_global==mmin(lerr_global)) = 0;
imagesc(sublist*100, rlist, lerr_global');
cm = gray(256); cm = cm(end:-1:1,:);
colormap(cm);
colorbar;
set(gca, 'FontSize', fs);
ylabel('Complexity \eta');
xlabel('K/N, % of kept eigenvectors');
if save_png
    saveas(gcf, [rep name '-all-error.png'], 'png');
end
if save_eps
    saveas(gcf, [rep name '-all-error.eps'], 'epsc');
end

