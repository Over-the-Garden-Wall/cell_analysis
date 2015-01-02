% close all
C = get_constants;
cell_nums = C.type.j(2);
cell_nums_minus1 = cell_nums(cell_nums ~= 15018);

all_mean_axis = [-0.6059    0.7955];
% all_mean_axis = [0.7955 -0.6059]
mean_axis_minus1 = [-0.6396    0.7687];

%     Mean (by angle):    -0.6059    0.7955
%     Mean sans 15018:    -0.6396    0.7687

quartile_data = cell(4,1);
bins = cell(4,1);



% [quartile_data{1} bins{1}] = get_stratification_by_dist(cell_nums);
% [quartile_data{2} bins{2}] = get_stratification_by_dist(cell_nums_minus1);
[quartile_data{3} bins{3}] = get_stratification_by_dist(cell_nums, all_mean_axis);
% [quartile_data{4} bins{4}] = get_stratification_by_dist(cell_nums_minus1, mean_axis_minus1);
% [quartile_data{5} bins{5}] = get_stratification_by_dist(cell_nums, []);


for n = 3:3;
    figure;    
    ar = quartile_data{n}(3:5,:)';
    
    lim = find(isnan(ar(:,1)),1,'first');
    if isempty(lim)
        lim = size(ar,1);
    else
        lim = lim-1;
    end
    
    hold on;
    
    area(bins{n}, ar(1:lim,3), 'FaceColor', [.75 .75 1], 'lineStyle', 'none');
    area(bins{n}, ar(1:lim,1), 'FaceColor', [1 1 1], 'lineStyle', 'none');
    plot(bins{n}, ar(1:lim,2), 'b')
    
    
end