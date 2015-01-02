close all

C = get_constants;

[quartile_data bins full_x full_data] = get_stratification_by_dist(C.type.j, [-0.6396    0.7687]);

quart_diff = quartile_data;
for n = 2:7; 
    quart_diff(n,:) = quartile_data(n,:)-quartile_data(n-1,:); 
end

quart_diff(5,:) = quart_diff(5,:)+quart_diff(4,:);
quart_diff(4,:) = [];


cmap = [1, .9, .75, .6, .75 .9 0]'*[1 1 0];
cmap(:,3) = 1;


figure; gca; hold all;
area_h = area(bins, quart_diff', 'LineStyle', 'none');

for n = 1:6
    set(area_h(n), 'FaceColor', cmap(n,:));
end





% for n = 1:6
%     area_h(n) = area(bins, quart_diff(n,:)','LineStyle', 'none', 'FaceColor', cmap(n,:));
% end


line_h = plot(bins, quartile_data(4,:), 'Color', cmap(end,:));
% colormap(cmap);

legend([area_h(2:4) line_h], {'Full arbor breadth', 'Inner 80 percentile', 'Inner 50 percentile', 'Median depth'});
ylabel('IPL depth (%)');
xlabel('Distance along soma-distal axis (microns)');

title('Stratification of all J cells by depth');