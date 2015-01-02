function soma_point = get_soma_loc(cell_id)



    C = get_constants;

    soma_fn = [C.soma_dir '/cell_' num2str(cell_id) '_soma.mat'];
    
    if exist(soma_fn, 'file')
        load(soma_fn)
    else
        disp(['soma point for ' num2str(cell_id) ' not found, calculating...']);
    
        fn = [C.point_dir '/cell_' num2str(cell_id) '_surface.mat'];
        load(fn);


        C = get_constants;
        depth = C.f(surface_points(:,1));

        %is the cell body near 0 or 100?


        d_hist = hist(depth, C.strat_x);
        d_hist = d_hist/sum(d_hist);

        [max_val strat_loc] = max(d_hist);

        poss_trunk_locs = d_hist < .01 & d_hist > 0;
        min_counts = sum(poss_trunk_locs(1:strat_loc));
        max_counts = sum(poss_trunk_locs(strat_loc:end));

        if min_counts > max_counts
            min_depth = min(depth);
            soma_points = surface_points(depth < (min_depth + C.soma_loc_threshold),:);
        else
            max_depth = max(depth);
            soma_points = surface_points(depth > (max_depth - C.soma_loc_threshold),:);
        end
        soma_point = mean(soma_points,1);
        save(soma_fn, 'soma_point');
        
    end
        
    
end