function plot_wave_solution(c,F,options)

% plot_wave_solution
%
%   plot_wave_solution(c,F,nb);
%
%   c is the speed
%   F(:,i) is the solution at time i.
%
%   Copyright (c) 2007 Laurent Demanet and Gabriel Peyre


ndims = 2;
if size(c,1)==1 || size(c,2)==1
    ndims = 1;
end

options.null = 0;
nb = getoptions(options, 'nplot', 6);
nt = size(F,ndims+1);
ind = round( linspace(1,nt,nb) );
t = getoptions(options, 't', linspace(0,1,nt));

savefile = getoptions(options, 'savefile', []);

n = size(F,1);
x = linspace(0,1,n);
lw = 1.5;
fs = 20;

if ndims==1

    clf;
    if isempty(savefile)
        subplot(nb+1,1,1);
    end
    h = plot(x, c); axis tight; 
    if isempty(savefile)
        title('Speed function.');
    end
    set(h, 'LineWidth', lw);
    if not(isempty(savefile))
        set(gca, 'FontSize', fs);
        % saveas(gcf, [savefile '-speed.png']);
        saveas(gcf, [savefile '-speed.eps'], 'epsc');
    end
    for i=1:nb
        str = ['t=' num2str(t(ind(i)), 2)];
        if isempty(savefile)
            subplot(nb+1,1,1+i);
        else
            clf;
        end
        h = plot(x,F(:,ind(i))); axis tight;
        axis([0 1 mmin(F) mmax(F)+eps]);
        set(h, 'LineWidth', lw);
        if isempty(savefile)
            title(str);
        end
        if not(isempty(savefile))
            set(gca, 'FontSize', fs);
            % saveas(gcf, [savefile '-' num2str(i) '.png']);
            saveas(gcf, [savefile '-' num2str(i) '.eps'], 'epsc');
        end
    end

else
    p = 2;
    if nb>10
        p = 3;
    end
    q = ceil(nb/p);
    clf;
    imageplot(c, 'speed', p, q,1);
    for i=1:nb-1
        str = ['t=' num2str(t(ind(i)), 2)];
        imageplot(F(:,:,ind(i)), str, p, q,i+1);
    end
end