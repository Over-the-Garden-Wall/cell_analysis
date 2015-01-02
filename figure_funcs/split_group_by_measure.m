function [g, g_val, g_violation] = split_group_by_measure(cell_nums, val, violation_threshold)


    hull = {zeros(0,2), zeros(0,2)};
    
    [val, val_order] = sort(val);
    cell_nums = cell_nums(val_order);    
    
    g = {[], []};
    g_val = {[], []};
    g_violation = {[], []};
    left_out_vals = [];
    
    left_out = [];
    
    phase = 1;
    
    while ~isempty(cell_nums)
        if phase == 1;
            cell_ind = 1;
        else
            cell_ind = length(cell_nums);
        end
        
        c_d = cell_data(cell_nums(cell_ind));
        cdh = [];
        [cdh(:,1), cdh(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
        
        h = polybool('intersection', cdh(:,1), cdh(:,2), hull{phase}(:,1), hull{phase}(:,2));
        if ~isempty(h)
            [h(:,1), h(:,2)] = polybool('intersection', cdh(:,1), cdh(:,2), hull{phase}(:,1), hull{phase}(:,2));
            intersect_area = poly_area(h);
        else
            intersect_area = 0;
        end
        
        if intersect_area / c_d.hull_area >= violation_threshold
            
            left_out(end+1) = cell_nums(cell_ind);
            left_out_vals(end+1) = val(cell_ind);
            
        else
        
            h = [];
            [h(:,1), h(:,2)] = polybool('union', cdh(:,1), cdh(:,2), hull{phase}(:,1), hull{phase}(:,2));
            hull{phase} = h;
            
            g{phase} = [g{phase} cell_nums(cell_ind)];
            g_violation{phase} = [g_violation{phase} intersect_area/c_d.hull_area];
            g_val{phase} = [g_val{phase} val(cell_ind)];
            
            phase = 3-phase;
        
        end
        
        cell_nums(cell_ind) = [];
        val(cell_ind) = [];
    end
    
    
    for k = 1:length(left_out)
        c_d = cell_data(left_out(k));
        cdh = [];
        [cdh(:,1), cdh(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
       
        h = [];
        
        intersect_area = zeros(2,1);
        for phase = 1:2
            h = polybool('intersection', cdh(:,1), cdh(:,2), hull{phase}(:,1), hull{phase}(:,2));
            if ~isempty(h)
                [h(:,1), h(:,2)] = polybool('intersection', cdh(:,1), cdh(:,2), hull{phase}(:,1), hull{phase}(:,2));
                intersect_area(phase) = poly_area(h);
            else
                intersect_area(phase) = 0;
            end
        end
    
        if intersect_area(1) < intersect_area(2)
            chosen_phase = 1;
        else
            chosen_phase = 2;
        end
        
        h = [];
        [h(:,1), h(:,2)] = polybool('union', cdh(:,1), cdh(:,2), hull{chosen_phase}(:,1), hull{chosen_phase}(:,2));
        hull{chosen_phase} = h;

        g{chosen_phase} = [g{chosen_phase} left_out(k)];
        g_violation{chosen_phase} = [g_violation{chosen_phase} intersect_area(chosen_phase)/c_d.hull_area];
        g_val{chosen_phase} = [g_val{chosen_phase} left_out_vals(k)];
    end
end