C = get_constants;
close all
cell_nums = unique([C.type.on_bc C.type.off_bc]);
% cell_nums = [];

% curr_types = {'t1', 't2', 't3a', 't3b', 't4'};
% curr_lbls = curr_types;

% for k = 1:length(curr_types);
%     cell_nums = [cell_nums C.type.(curr_types{k})];
% end
num_cells = length(cell_nums);


cutoff = 0;

depth_data = zeros(num_cells,3);

critical_depths = [.25, .5, .75, .85, .10];
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
depth_MSE = zeros(num_cells,1);
less_than_80 = zeros(num_cells,1);
num_skele_leaves = zeros(num_cells,1);
num_skele_branches = zeros(num_cells,1);
RF_displacement = zeros(num_cells,1);

early_peak = zeros(num_cells,1);
early_mean = zeros(num_cells,1);
early_min = zeros(num_cells,1);
trunk_width = zeros(num_cells,1);

curr_lbls = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5o', 'BC5i', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC', 'no id'};
curr_types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5o', 'BC5i', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC'};

% curr_lbls = {'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC', 'no id'};
% curr_types = {'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC'};


% [a b off_sac_conn_density] = connectivity_histogram(cell_nums, C.type.sure_off_sac, 1000);
for k = 1:num_cells;
    try
        c_d = cell_data(cell_nums(k));
        p = c_d.get_surface;
%         s = load([C.skele_dir 's' num2str(cell_nums(k)) '.mat']);
%         p = s.nodes;
        
        d = C.f(p(:,1));
        
        if sum(d>50) > length(d)*.25
            cutoff = 40;
        else
            cutoff = 0;
        end
        
        early_d = d(d>=10 & d<20,:);
        early_mean(k) = size(early_d,1);
        
        early_hist = hist(early_d, 10:.05:20);
        early_peak(k) = max(early_hist);
        early_min(k) = min(early_hist);
        
        load([C.skele_dir 's' num2str(cell_nums(k)) '.mat']);
        node_d = C.f(nodes(:,1));
        trunk_width(k) = mean(node_diameter(node_d > 0 & node_d < 10));
        
        
        midp = c_d.get_midpoint;
%         x = pca(p(d<cutoff,:));
%         trunk_direction(k,:) = x(:,1)';
        
        

        if isnan(cutoff)
            s = strat_notrunk(cell_nums(k));
            for n = find(s'==0);
                d(round(d)== n + C.strat_x(1) - 1) = [];
            end
            s;
            
        else
            p(d<cutoff,:) = [];
            d(d<cutoff) = [];
        end
        
        
        
%         
%         sixty_d = d(round(length(d)*.4));
%         p(d<=sixty_d,:) = [];
        hull_p = convhull(p(:,1), p(:,2));    
        hull = [];
        [hull(:,1), hull(:,2)] = poly2cw(p(hull_p,1), p(hull_p,2));
        
        hull_areas(k) = polyarea(hull(:,1), hull(:,2));
%         hull_areas(k) = c_d.hull_area;
        
        
        
        d = sort(d);
        d_len = length(d);
        
        for n = 1:length(critical_depths)
            depth_data(k,n) = d(round(d_len*critical_depths(n)));
        end
        
        [strat_peak(k) mode_ind] = max(c_d.stratification);
        mode_depth(k) = C.strat_x(mode_ind);
        total_points(k) = c_d.V * prod(C.res) / 1000^3;
        
        chat_split(k) = mean(d<64);
        less_than_80(k) = mean(d<80);
        strat_mean(k) = mean(d);
        depth_MSE(k) = mean((d - strat_mean(k)).^2);
        
        p = round(p(:,2:3)/1000);
        p = unique(p,'rows');
        coarse_area(k) = size(p,1);
        
        beyond_mode(k) = mean(d>mode_depth(k));
        
        for m = 1:length(curr_types)
            if any(cell_nums(k) == C.type.(curr_types{m}))
                curr_class(k) = m;
            end
        end            
        
        
%         num_nodes = size(nodes,1);
%         is_branch = false(num_nodes,1);
%         is_leaf = false(num_nodes,1);
% 
% 
%         for m = 1:size(nodes,1)
%             is_branch(m) = sum(edges(:)==m) >= 3;
%             is_leaf(m) = sum(edges(:)==m) == 1;
%         end
%         [dummy, top_node] = max(nodes(:,1));
%         
%         num_skele_leaves(k) = sum(is_leaf);
%         num_skele_branches(k) = sum(is_branch);
%         RF_displacement(k) = sum((midp(2:3)-nodes(top_node,(2:3))).*C.ventral_axis);
        
        
    catch ME
        disp(ME.message);
    end
    
%     if depth_data(k,3)-depth_data(k,1)>10
%         figure; plot(c_d.stratification); title(num2str(cell_nums(k)));
%     end
    


end;

is_valid = depth_data(:,1)~=0;

curr_class(~is_valid) = 0;

% depth_data = depth_data(is_valid,:);
% hull_areas = hull_areas(is_valid);
% cell_nums = cell_nums(is_valid);


% 
% figure; scatter(depth_data(:,1), depth_data(:,3));
% 
% hold on; plot([59 59], [min(depth_data(:,3)) max(depth_data(:,3))], 'r');
% hold on; plot([47 47], [min(depth_data(:,3)) max(depth_data(:,3))], 'r');
g = cell(14,1);
% 

%1-5 is off
g{1} = depth_data(:,3) < 50;
g{2} = depth_data(:,3) > 50;

g{3} = g{2} & depth_data(:,4) < 65;
g{4} = g{2} & depth_data(:,4) > 76.5;
g{2} = g{2} & ~g{3} & ~g{4};
g{5} = g{3} & depth_data(:,1) < 52.5;
g{3} = g{3} & ~g{5};
g{6} = g{5} & (depth_data(:,4) - depth_data(:,1)) < 9.75;
g{5} = g{5} & ~g{6};
g{7} = g{3} & coarse_area(:) > 325;
g{3} = g{3} & ~g{7};
g{8} = g{4} & coarse_area(:) < 90;
g{9} = g{4} & coarse_area(:) > 200;
g{4} = g{4} & ~g{8} & ~g{9};

g_labels = {'BC5t', 'BC5o', 'BC5i', 'XBC', 'BC6', 'BC7', 'BC89', 'RBC'};
g = {g{5}, g{6}, g{3}, g{7}, g{4}, g{2}, g{9}, g{8}};

M = zeros(8);
for n = 1:8
    my_nums = cell_nums(g{n});
   for m = 1:8
       type_nums = C.type.(g_labels{m});
       for c = type_nums;
           if any(c==my_nums)
               M(n,m) = M(n,m) + 1;
           end
       end
       
   end
end

disp(M);
% 
% g{2} = depth_data(:,3) < 66 & depth_data(:,3) >= 40; 
% g{1} = depth_data(:,3) > 66;
% 
% %g{3} is t5w
% g{3} = depth_data(:,3) - depth_data(:,1) > 7.2 & g{2};
% g{2} = g{2} & ~g{3};
% 
% %g{4} is xbc
% g{4} = g{2} & hull_areas > 8.5*10^8;
% g{2} = g{2} & ~g{4};
% 
% %
% % g{5} = g{2} & beyond_mode > .5;
% % g{2} = g{2} & ~g{5};
% %
% 
% %other types 5s
% g{5} = g{2} & depth_data(:,1) < mean(depth_data(g{2},1));
% g{2} = g{2} & ~g{5};
% 
% %g{6} is t7
% g{6} = g{1} & depth_data(:,3) < 75 & depth_data(:,3)-depth_data(:,1) < 10;
% g{1} = g{1} & ~g{6};
% 
% %g{7} is t89
% g{7} = g{1} & hull_areas > 5*10^8;
% g{1} = g{1} & ~g{7};
% 
% %g{8} is t6
% g{8} = g{1} & hull_areas > 1*10^8; %& less_than_80 > .4;
% g{1} = g{1} & ~g{8}; 
%g{1} is RBC


class_summary = [curr_class, zeros(length(curr_class),1)];
for k = 1:length(g);
    class_summary(g{k},2) = k;
end

M = zeros(max(class_summary(:,1))+1, max(class_summary(:,2))+1);
for n = 1:length(curr_class); 
    M(class_summary(n,1)+1, class_summary(n,2)+1) = M(class_summary(n,1)+1, class_summary(n,2)+1) + 1; 
end

% cmap = make_colormap(1-min(curr_class)+max(curr_class), 4);
% cmap(end, :) = [0 0 0];
% cmap(4, :) = [0 0 0];
% cmap(max(curr_class), :) = [0 .5 0];
cmap = [3 61 0; 215 148 0; 0 186 247; 166 217 125; 120 0 255; 56 255 18; 252 30 119; 212 38 255; 247 255 10; 0 179 84; 255 33 33; 255 71 3; 25 165 212; 150 100 250; 0 0 0]/255;


eval_strs = {'[depth_data(curr_class==k,4), coarse_area(curr_class==k)]', ...
    '[depth_data(curr_class==k,1), depth_data(curr_class==k,4)]', ...
    '[depth_data(curr_class==k,4), depth_data(curr_class==k,4)-depth_data(curr_class==k,1)]', ...
    '[depth_data(curr_class==k,1), depth_data(curr_class==k,4)-depth_data(curr_class==k,1)]', ...    
    '[depth_data(curr_class==k,5), coarse_area(curr_class==k)]', ...   
    '[depth_data(curr_class==k,1), coarse_area(curr_class==k)]'};

close all;

for m = 1:length(eval_strs)

figure; hold all

h = [];
k_list = [1:max(curr_class) 0];

BC_set = 1:length(k_list);
BC_set = [4 5];

for kc = BC_set
    k = k_list(kc);
    plot_data = eval(eval_strs{m});
    
    h(kc) = scatter(plot_data(:,1), plot_data(:,2), '*', 'markerEdgeColor', cmap(kc,:));
%     t_list = find(curr_class'==k);
%     for t = 1:length(t_list)
%         tc = t_list(t);
%         text(plot_data(t,1), plot_data(t,2), num2str(cell_nums(tc)-60000));
%     end
%     
    
%     h(k) = ellipse_pca([plot_data(:,1), plot_data(:,2)]);
%     set(h(k), 'lineWidth', 2, 'lineStyle', '--', 'Color', cmap(k+1,:));
end
prep_figure(gcf,gca, 'xlabel', 'xaxis', 'ylabel', 'yaxis', 'legend', curr_lbls(BC_set));
% xlabel('25th perc'); ylabel('75th perc');
% legend(h(2:end), curr_lbls([ 1:end-1]));

end




for k = 1:length(g_labels)
    fprintf([g_labels{k} ' = ']);
    fprintf('%d, ', cell_nums(g{k}));
    fprintf('\n');
end



k_list = [1:max(curr_class) 0];
eval_strs = {'depth_data(:,4)', 1:5, 6:14, '85th Percentile Depth', 'OFF', 'ON'; ...
    'depth_data(:,4)', 1:2, 3:5, '85th Percentile Depth', 'BC1,2', 'BC3a,3b,4'; ...
    'depth_data(:,4) - depth_data(:,1)', 1, 2, 'Stratification Width', 'BC1', 'BC2'; ...
    'depth_data(:,5)', 5, 3:4, '10th Percentile Depth', 'BC4', 'BC3a,3b'; ...
    'early_peak(:,1)', 5, 3:4, 'early peak', 'BC4', 'BC3a,3b'; ...
    'early_mean(:,1)', 5, 3:4, 'early mean', 'BC4', 'BC3a,3b'; ...
    'early_min(:,1)', 5, 3:4, 'early min', 'BC4', 'BC3a,3b'; ...
    'trunk_width(:,1)', 5, 3:4, 'trunk', 'BC4', 'BC3a,3b'; ...
    'total_points(:,1)', 4, 3, 'Volume (um)', 'BC3b', 'BC3a'; ...    
    'depth_data(:,4)', 6:9, 10:14, '85th Percentile Depth', 'BC5 & XBC', 'BC6,7,8,9 & RBC'; ...
    'depth_data(:,1)', [6 8], [7 9], '25th Percentile Depth', 'BC5t,5o', 'BC5i & XBC'; ...
    'depth_data(:,4) - depth_data(:,1)', 8, 6, 'Stratification Width', 'BC5o', 'BC5t'; ...
    'coarse_area(:,1)', 7, 9, 'Lateral Coverage (um)', 'BC5i', 'XBC'; ...
    'depth_data(:,4)', 11, [10 12:14], '85th Percentile Depth', 'BC7', 'BC6,8,9 & RBC'; ...
    'coarse_area(:,1)', [10 14], 12:13, 'Lateral Coverage (um)', 'BC6 & RBC', 'BC8,9'; ...
    'coarse_area(:,1)', 14, 10, 'Lateral Coverage (um)', 'RBC', 'BC6'; ...    
    'depth_data(:,4)', 12, 13, '85th Depth', 'BC8', 'BC9'; ... 
    'depth_data(:,1)', [6 8], 7, '25th Depth', 'BC5ot', 'BC5i'; ... 
    'depth_data(:,4)', 6, [7 8], '85th Depth', 'BCo', 'BC5it'; ... 
    
    };

is_me = [];
sub_data = [];    

for n = 1:size(eval_strs,1)
    figure;
    data = eval(eval_strs{n,1});
    for l = 1:2
        is_me{l} = false(size(curr_class,1),1);
        for k = 1:length(eval_strs{n,1+l})
            is_me{l} = is_me{l} | curr_class == eval_strs{n,1+l}(k);
        end
        sub_data{l} = data(is_me{l});
        
    end
    
    sub_data{2}(end+1:length(sub_data{1})) = NaN;
    sub_data{1}(end+1:length(sub_data{2})) = NaN;
    
    hist([sub_data{1} sub_data{2}], round((length(sub_data{1})+length(sub_data{2}))/4));
    stack_bars(gca);
    prep_figure(gcf,gca,'xlabel', eval_strs{n,4}, 'legend', eval_strs(n,5:6));
%     set(gca, 'box', 'off');
%     saveas(gcf, [C.image_dir, 'class_hist_', num2str(n), '.eps']);
end
    



close all


k_list = [1:max(curr_class) 0];
% eval_strs = {'depth_data(:,4)', 6:14, '85th Percentile Depth', 1, [65 76.5]; ...
%     'depth_data(:,4) - depth_data(:,1)', 1:2, 'Stratification Width', .5, 13; ...
%     'total_points(:,1)', 3:5, 'Volume (um)', 1, 48; ...    
%     'depth_data(:,1)', 6:9, '25th Percentile Depth', .5, 52.75; ...
%     'depth_data(:,4) - depth_data(:,1)', 6:7, 'Stratification Width', 1/3, 9.5; ...    
%     'coarse_area(:,1)', 8:9, 'Lateral Coverage (um)', 10, 300; ...
%     'coarse_area(:,1)', [10 12:14], 'Lateral Coverage (um)', 10, [90 200]; ...
%     };


eval_strs = {'depth_data(:,4)', 6:14, '85th Percentile Depth', 1, []; ...
    'depth_data(:,4) - depth_data(:,1)', 1:2, 'Stratification Width', .5, []; ...
    'total_points(:,1)', 3:5, 'Volume (um)', 1, []; ...    
    'depth_data(:,1)', 6:8, '25th Percentile Depth', .5, []; ...
    'depth_data(:,4) - depth_data(:,1)', 6:7, 'Stratification Width', 1/3, []; ...    
    'coarse_area(:,1)', 8:9, 'Lateral Coverage (um)', 10, []; ...
    'coarse_area(:,1)', [10 12:14], 'Lateral Coverage (um)', 10, []; ...
    'coarse_area(:,1)', [10 12:14], 'Lateral Coverage (um)', 10, []; ...
    };

is_me = [];
sub_data = [];    

clrs = colormap('Lines');

for n = 1:size(eval_strs,1)
    
    data = eval(eval_strs{n,1});
%     hist_data = nan(size(data,1), length(eval_strs{n,2}));
%     
    is_valid = false(size(data,1),1);
    for l = 1:length(eval_strs{n,2});
        is_valid = is_valid | curr_class==eval_strs{n,2}(l);        
%         my_class = eval_strs{n,2}(l);
%         hist_data(1:sum(curr_class==my_class), l) = data(curr_class==my_class);
        
    end
    data = data(is_valid);
        
    hist_data = nan(size(data,1), length(eval_strs{n,5}) + 1);
    splits = [-Inf eval_strs{n,5} Inf];
    for l = 1:length(splits)-1
        is_me = data > splits(l) & data <= splits(l+1);
        hist_data(1:sum(is_me),l) = data(is_me);
    end
    figure;
    hist(hist_data, floor(min(hist_data(:))):eval_strs{n,4}:ceil(max(hist_data(:))) );
    stack_bars(gca, clrs);%, cmap(eval_strs{n,2},:));
    set(gca,'YTick', 0:10:50)
    prep_figure(gcf,gca,'xlabel', eval_strs{n,3});
    set(gcf, 'position', [0 0 640 320]);
%     set(gca, 'box', 'off');
    save_fig([C.image_dir, 'class_hist_', num2str(n), '.eps'], gcf);
end
    

