% close all

C = get_constants;

cell_nums = C.type.sure_off_sac(3);
% cell_nums(cell_nums==15018) = [];
% cell_nums(cell_nums==17177) = [];

bins = 7.5:15:127.5;
[quartile_data, bins, full_x, full_data] = get_stratification_by_dist(cell_nums, [], bins, [.25 .5 .75], true, false,.1);

% M = [];
% for n = cell_nums
%     [blah, blah, x, dat] = get_stratification_by_dist(n, [], bins, [.25 .5 .75], true, false,1);
%     N = zeros(numel(dat),4);
%     N(:,1) = n;
%     num_x = length(x);
%     for b = 1:length(bins)
%         N((b-1)*num_x + (1:num_x),2) = b;
%         N((b-1)*num_x + (1:num_x),3) = x;
%         N((b-1)*num_x + (1:num_x),4) = dat(:,b);
%     end
%         
%     M = [M; N];
% end
% csvwrite('sac_strat_data.txt',M);



quart_diff = [quartile_data(1,:); quartile_data(3,:)-quartile_data(1,:)];

bin_cutoff = 5;

is_in = bins>=bin_cutoff;
quart_diff = quart_diff(:,is_in);
quartile_data = quartile_data(:,is_in);
bins = bins(is_in);

cmap = [1, .75, 0]'*[1 1 1] + (1-[1, .75, 0])'*C.colormap(6,:);
% cmap(:,3) = 1;


% figure; axis_h = subplot(2,2,1); hold all;
figure; hold on; ax = gca;
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


legend([area_h(2) line_h], {'25-75th percentiles', 'Median depth'});
ylabel('IPL depth (%)', 'FontSize', 20);
xlabel('Distance for soma (microns)', 'FontSize', 20);

set(gca,'XAxisLocation', 'bottom');
set(gca,'YDir', 'reverse');

% title('Stratification of all J cells by depth');



% 
% types = {'t1','t2','t3a','t3b','t4'};
% relevant_x = C.strat_x(C.strat_x >= relevant_portion(1) & C.strat_x <= relevant_portion(2));
% strats = zeros(length(relevant_x),length(types));
% 
% % strat_ax = subplot(2,2,2);
% hold on
% 
% 
% ax2 = axes('Position',get(ax,'Position'));
% 
% set(ax2, 'XAxisLocation','top');
% set(ax2, 'Color','none');
% set(ax2, 'YDir', 'reverse');
% 
% 
% 
% % set(ax2, 'XLim', [0 10]);
% 
% hold on
% 
% M = [];
% for k = 1:length(types);
%     cell_nums = C.type.(types{k});
%     for n = 1:length(cell_nums)
%         cell_dat = cell_data(cell_nums(n));
%         p = cell_dat.get_surface;
%         p = C.f(p(:,1));
%         p = p(p>=relevant_portion(1) & p<= relevant_portion(2));
%         hist_data = hist(p,relevant_x);
%         strats(:,k) = strats(:,k) + hist_data'/length(cell_nums);
%         
%         N = zeros(length(hist_data),4);
%         N(:,1) = cell_nums(n);
%         N(:,2) = k;
%         N(:,3) = relevant_x;
%         N(:,4) = hist_data;
%         
%         M = [M; N];
%     end
% %     strat_norm(k) = sum(strats(:,k));
% %     strats(:,k) = strats(:,k)/strat_norm(k);
%     plot(ax2, strats(:,k) * 410 / 1000 / 1000 * 100, relevant_x, 'LineWidth', 2, 'Color', C.colormap(k,:));
% %     plotyy(relevant_x, strats(:,k), relevant_x
% end
% 
% csvwrite('bip_strats_dat.txt',M);
% set(ax2, 'XLim', [0 1000]);
% set(ax2, 'XTick', 0:200:1000);
% set(ax2, 'YTick', []);
% 
% xlabel('BC surface voxels');
% 
% set(gca, 'FontSize', 20);
%     

% set(ax,'XTick', [-1 1]);


% legend([area_h(2) line_h], {'25-75th percentiles', 'Median depth'});
% legend(types);

% t_density = [2233 3212 1866 3254 3005];
% 
% is_valid = full_x >= relevant_x(1) & full_x <= relevant_x(end);
% full_data = full_data(is_valid,is_in);
% overlap = zeros(length(bins),length(types));
% 
% if size(full_data,1) > size(strats,1)
%     full_data = full_data(1:size(strats,1),:);
% elseif size(full_data,1) < size(strats,1)
%     strats = strats(1:size(full_data,1),:);
% end
% 
% for k = 1:length(types)
%     for n = 1:length(bins);
%         overlap(n,k) = sum(strats(:,k).*full_data(:,n))*strat_norm(k)*t_density(k);
%     end
% end
% 
% ol_ax = subplot(2,2,3);
% plot(bins, overlap, 'LineWidth', 2);
% legend(types);
% 
% set(ol_ax, 'YTick', []);
% title('expected total cell overlap by distance')
% xlabel('Distance along soma-distal axis (percent cell length)');
% 
% 
% subplot(2,2,4);
% 
% % ol_ratio = zeros(length(bins),length(types)-1);
% % for k = 1:length(types)-1;
% %     ol_ratio(:,k) = overlap(:,k+1)./overlap(:,1);
% % end
% 
% plot(bins, overlap(:,3)./overlap(:,1), 'LineWidth', 2);
% % legend(types(2:end));
% title('ratio of type 3a to type 1 bipolar cell overlap')
% xlabel('Distance along soma-distal axis (percent cell length)');
% 
