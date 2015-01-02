types = {'t1', 't2', 't3a', 't3b', 't4'};
num_types = length(types);

C = get_constants;

perc_dat = cell(num_types,1);
cell_nums = cell(num_types,1);
num_cells = zeros(num_types,1);

for k = 1:num_types;
    cell_nums{k} = C.type.(types{k});
    num_cells(k) = length(cell_nums{k});       
    perc_dat{k} = zeros(num_cells(k), 100);
    for ck = 1:num_cells(k);
        cell_dat = cell_data(cell_nums{k}(ck));
        p = cell_dat.get_surface;
        p = C.f(p(:,1));
        p = sort(p);
        num_points = length(p);
        for n = 1:100
            perc_dat{k}(ck,n) = p(round(num_points/100*n));
        end
    end
end
        

figure; hold all;
plot_points = [25 75];

for k = 1:num_types
    scatter(perc_dat{k}(:,plot_points(1)), perc_dat{k}(:,plot_points(2)), '*', 'MarkerEdgeColor', C.colormap(k,:));
end


search_range = 8:.01:12;
search_size = length(search_range);
E = zeros(search_size,1);
for k = 1:search_size;
    E(k) = sum(perc_dat{1}(:,plot_points(2)) - perc_dat{1}(:,plot_points(1)) > search_range(k));
    E(k) = E(k) + sum(perc_dat{2}(:,plot_points(2)) - perc_dat{2}(:,plot_points(1)) < search_range(k));    
end

is_best = E == min(E);
best_count = zeros(size(E));
best_count(1) = is_best(1);
for k = 2:search_size
    if is_best(k)
        best_count(k) = best_count(k-1)+1;
    end
end

best_best = max(best_count);
for k = search_size-1:-1:1
    if best_count(k) > 0 && best_count(k+1) > 0
        best_count(k) = best_count(k+1);
    end
end
is_best_best = best_best==best_count;


first_min_E = find(is_best_best,1,'first');
last_min_E = find(is_best_best,1,'last');

best_range = search_range(round(last_min_E/2 + first_min_E/2));

plot([8; 24], [8;24] + best_range, 'Color', [0 0 0], 'LineWidth', 2);
plot(gca, [10; 35], [35.5 35.5], 'Color', [0 0 0], 'LineWidth', 2);

prep_figure(gcf, gca, 'xlabel', '25^t^h percentile', 'ylabel', '75^t^h percentile', 'legend', {'BC1','BC2','BC3a','BC3b','BC4'});

