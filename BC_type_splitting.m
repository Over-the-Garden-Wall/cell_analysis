C = get_constants;

cell_nums = [C.type.on_bc C.type.off_bc];
% cell_nums = [];

% curr_types = {'t1', 't2', 't3a', 't3b', 't4'};
% curr_lbls = curr_types;

% for k = 1:length(curr_types);
%     cell_nums = [cell_nums C.type.(curr_types{k})];
% end
num_cells = length(cell_nums);


cutoff = 0;



% critical_depths = [.25, .5, .75, .90, .10];
critical_depths = .05:.05:.95;
depth_percentile = zeros(num_cells,length(critical_depths));
width = zeros(num_cells,floor(length(critical_depths)/2));

hull_areas = zeros(num_cells,1);
total_points = zeros(num_cells,1);
labels = zeros(num_cells,1);



curr_lbls = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC', 'no id'};
curr_types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC'};

% curr_lbls = {'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC', 'no id'};
% curr_types = {'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC8', 'BC9', 'RBC'};



for k = 1:num_cells;
    try
        c_d = cell_data(cell_nums(k));
        p = c_d.get_surface;
        d = C.f(p(:,1));
%         midp = c_d.get_midpoint;
%         x = pca(p(d<cutoff,:));
%         trunk_direction(k,:) = x(:,1)';
        
        to_delete = false(size(d,1),1);
        [s, total_points(k)] = strat_notrunk(cell_nums(k));
        for n = find(s'==0);
            to_delete(round(d)== n + C.strat_x(1) - 1) = true;
        end
        
        p = p(~to_delete, 2:3);
        d = d(~to_delete);
               
        d = sort(d);
        d_len = length(d);
        
        for n = 1:length(critical_depths)
            depth_percentile(k,n) = d(round(d_len*critical_depths(n)));
        end
        for n = 1:length(critical_depths)/2
            width(k,n) = depth_percentile(k,end-n+1) - depth_percentile(k,n);
        end
        
         
        hull_p = convhull(p(:,1), p(:,2));    
        hull = [];
        [hull(:,1), hull(:,2)] = poly2cw(p(hull_p,1), p(hull_p,2));
        
        hull_areas(k) = polyarea(hull(:,1), hull(:,2));
        
        
%         p = round(p(:,2:3)/1000);
%         p = unique(p,'rows');
%         coarse_area(k) = size(p,1);
        for n = 1:length(curr_types)
            if any(C.type.(curr_types{n})==cell_nums(k))
                labels(k) = n;
            end
        end


    catch ME
        disp(ME.message);
    end
    

end;


features = [depth_percentile width hull_areas total_points];
num_features = size(features,2);

features = features(labels > 0, :);
labels(labels==0) = [];
num_examples = length(labels);

%expand features into binary set
pairs = [];
[pairs(:,1), pairs(:,2)] = ind2sub(max(labels)*[1 1], 1:max(labels)^2);
pairs(pairs(:,2) <= pairs(:,1),:) = [];
num_pairs = size(pairs,1);

separability = zeros([num_features, max(pairs)]);

for f = 1:num_features
    for p = 1:num_pairs
        is_relevant = labels == pairs(p,1) | labels == pairs(p,2);
        sub_f = features(is_relevant, f);
        sub_lbl = labels(is_relevant);
        
        [sub_f, sort_ord] = sort(sub_f);
        sub_lbl = sub_lbl(sort_ord) == sub_lbl(1);
        num_subsamps = length(sub_lbl);
        
        cum_lbl1 = cumsum(sub_lbl);
        cum_lbl1 = cum_lbl1/cum_lbl1(end);
        cum_lbl2 = cumsum(1-sub_lbl);
        cum_lbl2 = cum_lbl2/cum_lbl2(end);
                
        errors = min([1 - cum_lbl1' + cum_lbl2'; 1 + cum_lbl1' - cum_lbl2']);
        [fewest_errors split_point] = min(errors);
        
        separability(f, pairs(p,1), pairs(p,2)) = fewest_errors;
        separability(f, pairs(p,2), pairs(p,1)) = fewest_errors;
    end
end