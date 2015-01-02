close all

C = get_constants;

bins = 10:10:200;

[quartile_data bins full_x full_data] = get_stratification_by_dist(C.type.j, [], bins, [.25 .5 .75]);

quart_diff = [quartile_data(1,:); quartile_data(3,:)-quartile_data(1,:)];

bin_cutoff = 5;

is_in = bins>=bin_cutoff;
quart_diff = quart_diff(:,is_in);
quartile_data = quartile_data(:,is_in);
bins = bins(is_in);

cmap = [1, .75, 0]'*[1 1 0];
cmap(:,3) = 1;


figure; axis_h = subplot(1,2,1); hold all;
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
xlabel('Distance along soma-distal axis (microns)');

title('Stratification of all J cells by depth');




types = {'t1','t2','t3a','t3b','t4'};
relevant_x = C.strat_x(C.strat_x >= relevant_portion(1) & C.strat_x <= relevant_portion(2));
strats = zeros(length(relevant_x),length(types));

strat_ax = subplot(1,2,2);
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
    strats(:,k) = strats(:,k)/sum(strats(:,k));
    plot(strats(:,k), relevant_x);
end

set(strat_ax,'XTick', [-1 1]);
legend(types);
