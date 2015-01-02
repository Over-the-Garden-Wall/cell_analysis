type = 't2',
c = colormap('Lines');
close all
C = get_constants;

cell_nums = C.type.(type);
num_cells = length(cell_nums);

hulls = cell(num_cells,1);
for k = 1:num_cells;
    cell_dat = cell_data(cell_nums(k));
    [hulls{k}(:,1), hulls{k}(:,2)] = poly2cw(cell_dat.hull_2d(:,1), cell_dat.hull_2d(:,2));
end

ol_hulls = cell(num_cells);

for k = 1:num_cells-1;
    for j = k+1:num_cells;
        t = polybool('intersection', hulls{k}(:,1), hulls{k}(:,2), hulls{j}(:,1), hulls{j}(:,2));
        if ~isempty(t)
            [ol_hulls{k,j}(:,1), ol_hulls{k,j}(:,2)] = polybool('intersection', hulls{k}(:,1), hulls{k}(:,2), hulls{j}(:,1), hulls{j}(:,2));
        end
    end
end

figure; hold all;


for k = 1:num_cells;
    plot(hulls{k}(:,1), hulls{k}(:,2), 'LineWidth', 2);
end

set(gca, 'XLim', [.4 2.4]*10^5);
set(gca, 'YLim', [.4 2.4]*10^5);
prep_figure(gcf,gca);



figure; hold all;
for k = 1:num_cells;   
    plot_cells(cell_nums(k), 1, .01, c(mod(k-1,64)+1,:));        
end

set(gca, 'XLim', [.4 2.4]*10^5);
set(gca, 'YLim', [.4 2.4]*10^5);
prep_figure(gcf,gca);