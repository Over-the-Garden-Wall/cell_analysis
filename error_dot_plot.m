function h = error_dot_plot(Y, Yerr, xlabels, clrs)

    
    if ~exist('xlabels', 'var') || isempty(xlabels);
        xlabels = cell(size(Y,1),1);
        for k = 1:length(xlabels)
            xlabels{k} = ['cat ' num2str(k)];
        end
    end
    
    if ~exist('clrs', 'var') || isempty(clrs);
        clrs = colormap('Lines');
    end
    
    h = zeros(size(Y,2),1);
    
    hold all;
    
    for n = 1:size(Y,2)
        h(n) = errorbar(1:size(Y,1), Y(:,n)', Yerr(:,n)', 'lineWidth', 2, 'color', clrs(n,:), 'marker', 'o', 'markerSize', 10, 'lineStyle', 'none');
        
    end
    set(gca, 'XTick', 1:size(Y,1), 'XTickLabel', xlabels, 'XLim', [0 size(Y,1)+1]);
end