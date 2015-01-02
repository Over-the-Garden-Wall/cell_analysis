close all

C = get_constants;
figure;

axis_lims = [-40 100 0 50];

for k = 1:length(C.type.j)
cell_num = C.type.j(k);
[quartile_data bins full_x full_data] = get_stratification_by_dist(cell_num);

quart_diff = quartile_data;
for n = 2:7; 
    quart_diff(n,:) = quartile_data(n,:)-quartile_data(n-1,:); 
end

quart_diff(5,:) = quart_diff(5,:)+quart_diff(4,:);
quart_diff(4,:) = [];


cmap = [1, .9, .75, .6, .75 .9 0]'*[1 1 0];
cmap(:,3) = 1;


subplot(7, 1,k); hold all;
area_h = area(bins, quart_diff', 'LineStyle', 'none');

for n = 1:6
    set(area_h(n), 'FaceColor', cmap(n,:));
end





% for n = 1:6
%     area_h(n) = area(bins, quart_diff(n,:)','LineStyle', 'none', 'FaceColor', cmap(n,:));
% end


line_h = plot(bins, quartile_data(4,:), 'Color', cmap(end,:));
axis(axis_lims);
% colormap(cmap);

% ylabel('IPL depth (%)');
% xlabel('Distance along soma-distal axis (microns)');

title([num2str(cell_num)]);
end
legend([area_h(2:4) line_h], {'Full arbor breadth', 'Inner 80 percentile', 'Inner 50 percentile', 'Median depth'});
