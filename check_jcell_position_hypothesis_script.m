function check_jcell_position_hypothesis_script
    
    C = get_constants;
    
    s_cell_nums = C.type.off_sac;
    num_s = length(s_cell_nums);
    j_cell_nums = C.type.j;
    num_j = length(j_cell_nums);
    legend_str = cell(num_j,1);
    
    s_pos = zeros(num_s,2);
    j_pos = zeros(num_j,2);
    
    for s = 1:num_s
        c = cell_data(s_cell_nums(s));
        temp_loc = c.get_midpoint(true);
        s_pos(s,:) = temp_loc(2:3);
    end
    for j = 1:num_j
        c = cell_data(j_cell_nums(j));
        temp_loc = c.get_midpoint(true);
        j_pos(j,:) = temp_loc(2:3);
        legend_str{j} = num2str(c.cell_id);
    end
    
    
    scatter_handles = zeros(num_j,1);
    figure; hold all
    
    for j = 1:num_j;
        rel_s_pos = s_pos;
        for d = 1:2
            rel_s_pos(:,d) = rel_s_pos(:,d) - j_pos(j,d);
        end
%         figure; 
        scatter_handles(j) = scatter(rel_s_pos(:,1), rel_s_pos(:,2), '*');
        
%         title(legend_str{j});
    end
    
%     legend(scatter_handles, legend_str);
   
    figure;
    hold all;
    scatter(s_pos(:,1), s_pos(:,2), '*');
    scatter(j_pos(:,1), j_pos(:,2), 'o');
    

end