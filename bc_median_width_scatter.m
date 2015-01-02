C = get_constants;

cell_nums = C.type.on_bc;
% cell_nums = [];

% curr_types = {'t1', 't2', 't3a', 't3b', 't4'};
% curr_lbls = curr_types;

% for k = 1:length(curr_types);
%     cell_nums = [cell_nums C.type.(curr_types{k})];
% end
num_cells = length(cell_nums);

cutoff = 0;

depth_data = zeros(num_cells,3);

critical_depths = [.25, .5, .75, .90];
hull_areas = zeros(num_cells,1);
mode_depth = zeros(num_cells,1);
total_points = zeros(num_cells,1);
chat_split = zeros(num_cells,1);
coarse_area = zeros(num_cells,1);
beyond_mode = zeros(num_cells,1);
trunk_direction = zeros(num_cells,3);
strat_peak = zeros(num_cells,1);
strat_mean = zeros(num_cells,1);
curr_class = zeros(num_cells,1);


curr_types = {'t5w', 't5l', 't5h', 'xbc', 't6', 't7', 't89', 'tRBC'};
curr_lbls = {'BC5w', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC89', 'RBC'};




for k = 1:num_cells;
    try
        c_d = cell_data(cell_nums(k));
        p = c_d.get_surface;
        d = C.f(p(:,1));
        
%         x = pca(p(d<cutoff,:));
%         trunk_direction(k,:) = x(:,1)';
        
        d(d<cutoff) = [];
        d = sort(d);
        d_len = length(d);
        
        for n = 1:length(critical_depths)
            depth_data(k,n) = d(round(d_len*critical_depths(n)));
        end
        
        hull_areas(k) = c_d.hull_area;
        [strat_peak(k) mode_ind] = max(c_d.stratification);
        mode_depth(k) = C.strat_x(mode_ind);
        total_points(k) = c_d.V;
        
        chat_split(k) = mean(d<64);
        
        strat_mean(k) = mean(d);
        
        p = round(p(:,2:3)/1000);
        p = unique(p,'rows');
        coarse_area(k) = size(p,1);
        
        beyond_mode(k) = mean(d>mode_depth(k));
        
        for m = 1:length(curr_types)
            if any(cell_nums(k) == C.type.(curr_types{m}))
                curr_class(k) = m;
            end
        end
        
    catch ME
        disp(ME.message);
    end
    
%     if depth_data(k,3)-depth_data(k,1)>10
%         figure; plot(c_d.stratification); title(num2str(cell_nums(k)));
%     end
    
end;

is_valid = depth_data(:,1)~=0;

depth_data = depth_data(is_valid,:);
hull_areas = hull_areas(is_valid);
cell_nums = cell_nums(is_valid);
% 
% figure; scatter(depth_data(:,1), depth_data(:,3));
% 
% hold on; plot([59 59], [min(depth_data(:,3)) max(depth_data(:,3))], 'r');
% hold on; plot([47 47], [min(depth_data(:,3)) max(depth_data(:,3))], 'r');

% 
% g{1} = depth_data(:,1) < 47;
% sg{1} = depth_data(:,1) > 59;
% 
% 
% sg{2} = ~sg{1} & ~g{1};
% 
% 
% g{9} = sg{2} & depth_data(:,3) > 70;
% sg{3} = sg{2} & ~g{9};
% 
% % fprintf('too shallow group: '); fprintf('%d ', cell_nums(too_shallow_g));
% % fprintf('\nshallow group: '); fprintf('%d ', cell_nums(shallow_g));
% % fprintf('\ndeep group: '); fprintf('%d ', cell_nums(deep_g));
% 
% 
% 
% g{2} = sg{3} & hull_areas > 8*10^8;
% sg{3} = sg{3} & ~g{2};
% 
% g{3} = sg{3} & ((depth_data(:,3)-depth_data(:,1)) > 7.5);
% sg{4} = sg{3} & ~g{3};
% g{4} = sg{4} & depth_data(:,1) > 54;
% g{5} = sg{4} & ~g{4};
% g{4} = g{4} & ~g{2};
% 
% 
% g{6} = sg{1} & depth_data(:,3) < 73;
% 
% 
% sg{5} = sg{1} & ~g{6};
% 
% 
% g{7} = sg{5} & hull_areas > 5*10^8;
% 
% sg{6} = sg{5} & ~g{7};
% 
% figure; scatter(depth_data(sg{6},2), total_points(sg{6}));
% 
% % g{8} = depth_data(:,2) > 85 & sg{6};
% g{8} = depth_data(:,2) > 77 & sg{6};
% % g{9} = total_points > 8.5*10^5 & sg{6} & ~g{8};
% % g{10} = sg{6} & ~g{8} & ~g{9};
% g{9} = (sg{6} & ~g{8}) | g{9};
% for k = 1:length(g);
%     g_name{k} = ['group ' num2str(k)];
% end
% 
% figure; hold all
% for k = 1:length(g)
%     scatter(depth_data(g{k},1), depth_data(g{k},3), '*');
% end
% xlabel('25th perc'); ylabel('75th perc');
% legend(g_name);
% 
% 
% figure; hold all
% for k = 1:length(g)
%     scatter(depth_data(g{k},4)-depth_data(g{k},2), hull_areas(g{k}), '*');
% end
% xlabel('stratification width (90-50)'); ylabel('hull area');
% legend(g_name);
% 
% figure; hold all
% for k = 1:length(g)
%     scatter(depth_data(g{k},2), total_points(g{k}), '*');
% end
% xlabel('50th perc (median)'); ylabel('cell volume');
% legend(g_name);
% 
% for k = 1:length(g);
%     fprintf('group %d: ', k); 
%     fprintf('%d ', cell_nums(g{k}));
%     fprintf('\n');
% end
% 

figure; hold all
for k = [5 7 8]
    is_me = curr_class==k;
    scatter(depth_data(is_me,3), hull_areas(is_me), '*');
end
legend(curr_lbls([5 7 8]));
xlabel('75th percentile depth'); ylabel('coverage');

% 
% figure; hold all
% for k = 1:4
%     is_me = curr_class==k;
%     scatter(depth_data(is_me,1), hull_areas(is_me), '*');
% end
% legend(curr_lbls);
% xlabel('25th percentile depth'); ylabel('Coverage');
% 
% 
% figure; hold all
% for k = 1:3
%     is_me = curr_class==k;
%     scatter(depth_data(is_me,1), depth_data(is_me,3)-depth_data(is_me,1), '*');
% end
% legend(curr_lbls);
% xlabel('25th percentile depth'); ylabel('75th - 25th Percentile');


% 
% figure; hold all
% for k = 4:5
%     scatter(depth_data(g{k},1), mode_depth(g{k}), '*');
% end
% 
% for k = 4:5
%     ck = find(g{k});
%     
%     for n = 1:length(ck)
%         text(depth_data(ck(n),1), mode_depth(ck(n)), num2str(cell_nums(ck(n))-60000));
%     end
% end
% 
% legend(g_name(4:5));
% 
% 
% figure; hold all
% for k = 3:5
%     scatter(depth_data(g{k},1), strat_peak(g{k}), '*');
% end
% legend(g_name(3:5));
% 
% 
% % for k = 4:5
% %     ck = find(g{k});
% %     
% %     for n = 1:length(ck)
% %         text(depth_data(ck(n),1), mode_depth(ck(n)), num2str(cell_nums(ck(n))-60000));
% %     end
% % end
% % xlabel('25th percentile depth'); ylabel('mode depth');
% % legend(g_name(4:5));
% % 
% % 
% % figure; hold all
% % for k = 7:9
% %     scatter(depth_data(g{k},2), hull_areas(g{k}), '*');
% % end
% % xlabel('50th percentile depth'); ylabel('hull area');
% % legend(g_name(7:9));
% % 
% % 
% % figure; hold all
% % for k = 7:9
% %     scatter(depth_data(g{k},4), total_points(g{k}), '*');
% % end
% % xlabel('90th percentile depth'); ylabel('volume');
% % legend(g_name(7:9));
% % 
% % 
% 
% 
% 
% % close all
% 
% figure; hist(depth_data(:,1), 20);
% prep_figure(gcf,gca,'xlabel', '25th percentile', 'title', 't5, xbc and misfit t6 vs t6+ split');
% 
% figure; hist(depth_data(sg{2},3), 20);
% prep_figure(gcf,gca,'xlabel', '75th percentile', 'title', 't5, xbc vs misfit t6 split');
% 
% figure; hist(coarse_area(sg{2}), 20);
% prep_figure(gcf,gca,'xlabel', 'hull area', 'title', 't5 vs xbc split');
% 
% figure; hist((depth_data(g{3} | g{4} | g{5},3)-depth_data(g{3} | g{4} | g{5},1)), 20);
% prep_figure(gcf,gca,'xlabel', '75th percentile - 25th percentile', 'title', 't5lh vs t5w split');
% 
% figure; hist(depth_data(g{4} | g{5},1), 20);
% prep_figure(gcf,gca,'xlabel', '25th percentile', 'title', 't5a vs t5b split?');
% 
% figure; hist(depth_data(sg{1},3), 20);
% prep_figure(gcf,gca,'xlabel', '75th percentile - 25th percentile', 'title', 't7 vs other t6+ split');
% 
% figure; hist(coarse_area(sg{1} & ~g{6}), 20);
% prep_figure(gcf,gca,'xlabel', 'hull area', 'title', 't8+9 vs RBC/t6');
% 
% figure; hist(depth_data(sg{6}, 2), 20);
% prep_figure(gcf,gca,'xlabel', 'median depth', 'title', 't6 vs RBC?');
