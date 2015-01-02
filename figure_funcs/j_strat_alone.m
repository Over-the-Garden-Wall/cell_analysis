% close all

C = get_constants;
% types = {'t1','t2','t3a','t3b','t4', 'sure_off_sac'};

cell_nums = C.type.j;
% bins = C.j_bins;
% bins = 20:20:180;
bins = 20:20:200;

% cell_nums(cell_nums==15018) = [];
% cell_nums(cell_nums==17177) = [];

[quartile_data bins full_x full_data] = get_stratification_by_dist(cell_nums, [], bins, [.25 .5 .75], true, false,.1);

quart_diff = [quartile_data(1,:); quartile_data(3,:)-quartile_data(1,:)];

bin_cutoff = 10;

is_in = bins>=bin_cutoff;
quart_diff = quart_diff(:,is_in);
quartile_data = quartile_data(:,is_in);
bins = bins(is_in);

c = colormap('Lines');

cmap = [1, .45, 0]'*[1 1 1] + (1-[1, .45, 0])'*c(1,:);
% cmap(:,3) = 1;


% figure; axis_h = subplot(2,2,1); hold all;
% figure; hold on; ax = gca;
area_h = area(ax, bins, quart_diff', 'LineStyle', 'none');

for n = 1:2
    set(area_h(n), 'FaceColor', cmap(n,:));
end


relevant_portion = [0 60];

set(ax,'YLim', relevant_portion);



% for n = 1:6
%     area_h(n) = area(bins, quart_diff(n,:)','LineStyle', 'none', 'FaceColor', cmap(n,:));
% end


line_h = plot(bins, quartile_data(2,:), 'LineWidth', 2, 'Color', cmap(end,:));
% colormap(cmap);

    set(gcf, 'Position', [0 0 640 640]);
    set(gca, 'Position', [.1 .1 .8 .8]);

set(gca, 'FontSize', 20);
% set(gca, 'YDir', 'reverse');
set(gca, 'XLim', [0 150]);
set(gca, 'XTick', 0:50:1000);

legend([area_h(2) line_h], {'25-75th percentiles', 'Median depth'});
ylabel('IPL depth (%)', 'FontSize', 20);
xlabel('Distance from soma (microns)', 'FontSize', 20);

set(gca,'XAxisLocation', 'bottom');
% set(gca,'YDir', 'reverse');

% title('Stratification of all J cells by depth');


