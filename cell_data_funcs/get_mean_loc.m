function mean_point = get_mean_loc(cell_id)



    C = get_constants;

    mean_fn = [C.soma_dir '/cell_' num2str(cell_id) '_mean.mat'];
    
    if exist(mean_fn, 'file')
        load(mean_fn)
    else
        disp(['mean point for ' num2str(cell_id) ' not found, calculating...']);
    
        fn = [C.point_dir '/cell_' num2str(cell_id) '_surface.mat'];
        load(fn);

        mean_point = mean(surface_points,1);
        save(mean_fn, 'mean_point');
        
    end
        
    
end