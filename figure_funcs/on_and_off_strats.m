function on_and_off_strats

    % close all

C = get_constants;

on_cells = C.type.on_sac;
off_cells = C.type.off_sac;

bins = C.sac_bins;
bins = [0:20:160];


[on_quartile_data, bins] = get_stratification_by_dist(on_cells, [], bins, [.25 .5 .75], true, false,.1);
[off_quartile_data, bins] = get_stratification_by_dist(off_cells, [], bins, [.25 .5 .75], true, false,.1);


is_in = ~isnan(on_quartile_data(1,:)) & ~isnan(off_quartile_data(1,:));
bins = bins(is_in);
on_quartile_data = on_quartile_data(:,is_in);
off_quartile_data = off_quartile_data(:,is_in);

on_quart_diff = [on_quartile_data(1,:); on_quartile_data(3,:)-on_quartile_data(1,:)];
off_quart_diff = [off_quartile_data(1,:); off_quartile_data(3,:)-off_quartile_data(1,:)];

on_c = [1 0 0];
off_c = [0 0 1];

on_area_c = [1 .5 .5];
off_area_c = [.5 .5 1];


% figure; axis_h = subplot(2,2,1); hold all;
figure; hold on; ax = gca;
on_area_h = area(ax, bins, on_quart_diff', 'LineStyle', 'none');
off_area_h = area(ax, bins, off_quart_diff', 'LineStyle', 'none');


set(on_area_h(1), 'FaceColor', [1 1 1]);
set(off_area_h(1), 'FaceColor', [1 1 1]);
set(on_area_h(2), 'FaceColor', on_area_c);
set(off_area_h(2), 'FaceColor', off_area_c);


relevant_portion = [0 100];

set(ax,'YLim', relevant_portion);



% for n = 1:6
%     area_h(n) = area(bins, quart_diff(n,:)','LineStyle', 'none', 'FaceColor', cmap(n,:));
% end

off_line_h = plot(bins, off_quartile_data(2,:), 'LineWidth', 2, 'Color', off_c);
on_line_h = plot(bins, on_quartile_data(2,:), 'LineWidth', 2, 'Color', on_c);
% colormap(cmap);

    set(gcf, 'Position', [0 0 640 640]);
    set(gca, 'Position', [.1 .1 .8 .8]);

set(gca, 'FontSize', 20);


legend([on_line_h on_area_h(2) off_line_h off_area_h(2)], {'ON-Median', 'ON-25 to 75th %', 'OFF-Median', 'OFF-25 to 75th %'});
ylabel('IPL depth (%)', 'FontSize', 20);
xlabel('Distance for soma (microns)', 'FontSize', 20);

set(gca,'XAxisLocation', 'bottom');
set(gca,'YDir', 'reverse');

% title('Stratification of all J cells by depth');


