function overlap_by_dist = contact_prob_depthXdist(target_cells, contact_cells, max_dist)
    

    target_dist = count_depth_by_distance(target_cells, max_dist);
    joint_dist = zeros(size(target_dist));
    
    [norm_strat total_strat] = get_mean_stratification(contact_cells, 1, 100);
    
    for n = 1:max_dist
        joint_dist(:,n) = target_dist(:,n).*total_strat(:);
    end
    
    overlap_by_dist = sum(joint_dist);
            
end