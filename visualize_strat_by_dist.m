
type = 'minij';

cell_nums = C.type.(type);
bins = C.([type '_bins']);

x_range = 0:100;

num_cells = length(cell_nums);

for k = 1:num_cells
    figure;

    [quartile_data bins full_x full_data] = get_stratification_by_dist(cell_nums(k), [], bins, [.25 .5 .75], true, false, 1);

    if full_x(1) < x_range(1);
        num_missed = x_range(1)-full_x(1);
        full_data = [zeros(num_missed,size(full_data,2)); full_data];
    elseif full_x(1) > x_range(1);
        num_missed = full_x(1) - x_range(1);
        full_data(1:num_missed,:) = [];
    end
    
    
    if full_x(end) < x_range(end);
        num_missed = x_range(end)-full_x(end);
        full_data = [full_data; zeros(num_missed,size(full_data,2))];
    elseif full_x(end) > x_range(end);
        num_missed = full_x(end) - x_range(end);
        full_data = full_data(1:end-num_missed,:);
    end
    
    vox_per_bin = zeros(1,size(full_data,2));
    for n = 1:size(full_data,2)
        vox_per_bin(n) = sum(full_data(:,n));
        full_data(:,n) = full_data(:,n)/vox_per_bin(n);
    end
    
        
        
    subplot(2,1,1); imagesc(full_data); title(num2str(cell_nums(k)));
    subplot(2,1,2); plot(vox_per_bin, 'LineWidth', 2);
end