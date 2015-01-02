close all

C = get_constants;

cell_nums = C.type.j;
% cell_nums(cell_nums==15018) = [];
cell_nums(cell_nums==17177) = [];


[quartile_data bins full_x full_data] = get_stratification_by_dist(cell_nums, [], 10:10:90, [.25 .5 .75], true, true);

quart_diff = [quartile_data(1,:); quartile_data(3,:)-quartile_data(1,:)];

bin_cutoff = 5;

is_in = bins>=bin_cutoff;
quart_diff = quart_diff(:,is_in);
quartile_data = quartile_data(:,is_in);
bins = bins(is_in);

cmap = [1, .75, 0]'*[1 1 0];
cmap(:,3) = 1;


figure; axis_h = subplot(2,2,1); hold all;
area_h = area(bins, quart_diff', 'LineStyle', 'none');

for n = 1:2
    set(area_h(n), 'FaceColor', cmap(n,:));
end


relevant_portion = [0 60];

set(axis_h,'YLim', relevant_portion);



% for n = 1:6
%     area_h(n) = area(bins, quart_diff(n,:)','LineStyle', 'none', 'FaceColor', cmap(n,:));
% end


line_h = plot(bins, quartile_data(2,:), 'Color', cmap(end,:));
% colormap(cmap);

legend([area_h(2) line_h], {'25-75th percentiles', 'Median depth'});
ylabel('IPL depth (%)');
xlabel('Distance along soma-distal axis (percent cell length)');

title('Stratification of all J cells by depth');




types = {'t1','t2','t3a','t3b','t4', 'off_sac'};
relevant_x = C.strat_x(C.strat_x >= relevant_portion(1) & C.strat_x <= relevant_portion(2));
strats = zeros(length(relevant_x),length(types));

strat_ax = subplot(2,2,2);
hold all



for k = 1:length(types);
    cell_nums = C.type.(types{k});
    for n = 1:length(cell_nums)
        cell_dat = cell_data(cell_nums(n));
        p = cell_dat.get_surface;
        p = C.f(p(:,1));
        p = p(p>=relevant_portion(1) & p<= relevant_portion(2));
        strats(:,k) = strats(:,k) + hist(p,relevant_x)';
    end
    strat_norm(k) = sum(strats(:,k));
    strats(:,k) = strats(:,k)/strat_norm(k);
%     plot(strats(:,k), relevant_x, 'LineWidth', 2);
end

grouped_strats = [mean(strats(:,1:2),2) mean(strats(:,3:5),2) strats(:,6)];
for k = 1:size(grouped_strats,2)
    plot(grouped_strats(:,k), relevant_x, 'LineWidth', 2);
end

set(strat_ax,'XTick', [-1 1]);
legend('types 1&2', 'types 3&4', 'off SACs');

t_density = [2233 3212 1866 3254 3005 1100];

is_valid = full_x >= relevant_x(1) & full_x <= relevant_x(end);
full_data = full_data(is_valid,is_in);
overlap = zeros(length(bins),length(types));

for k = 1:length(types)
    for n = 1:length(bins);
        overlap(n,k) = sum(strats(:,k).*full_data(:,n))*strat_norm(k)*t_density(k)/length(C.type.(types{k}));
    end
end

grouped_overlap = [mean(overlap(:,1:2),2) mean(overlap(:,3:5),2) overlap(:,6)];

ol_ax = subplot(2,2,3);
plot(bins, grouped_overlap, 'LineWidth', 2);
legend('types 1&2', 'types 3&4', 'off SACs');

set(ol_ax, 'YTick', []);
title('expected total cell overlap by distance')
xlabel('Distance along soma-distal axis (percent cell length)');


subplot(2,2,4);

% ol_ratio = zeros(length(bins),length(types)-1);
% for k = 1:length(types)-1;
%     ol_ratio(:,k) = overlap(:,k+1)./overlap(:,1);
% end

ratios = [grouped_overlap(:,2)./grouped_overlap(:,1), ...
    grouped_overlap(:,1)./grouped_overlap(:,3), ...
    grouped_overlap(:,2)./grouped_overlap(:,3)];

plot(bins, ratios, 'LineWidth', 2);
% legend(types(2:end));
title('Overlap Ratios')
xlabel('Distance along soma-distal axis (percent cell length)');
legend('t3&4:t1&2','t1&2:sac','t3&4:sac');
