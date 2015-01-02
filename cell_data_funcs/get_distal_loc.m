function distal_point = get_distal_loc(cell_id, varargin)

    p = inputParser;    
    p.addRequired('cell_id', @isnumeric);
    p.addOptional('use_soma', check_to_use_soma(cell_id), @islogical);
    
    p.parse(cell_id, varargin{:});    
    s = p.Results;

    C = get_constants;

    dist_fn = [C.soma_dir '/cell_' num2str(s.cell_id) '_dist.mat'];
    if exist(dist_fn, 'file')
        load(dist_fn)
    else
    
    
    
    
        mean_point = get_mean_point(s.cell_id, s.use_soma);
        
        disp(['distal point for ' num2str(s.cell_id) ' not found, calculating...']);
    
        fn = [C.point_dir '/cell_' num2str(s.cell_id) '_surface.mat'];
        load(fn);

        num_points = size(surface_points,1);
        dist = sqrt(sum((surface_points(:,2:3) - ones(num_points,1)*mean_point(2:3)).^2,2));
        
        [max_dist, max_ind] = max(dist);
        distal_point = surface_points(max_ind,:);
        
        save(dist_fn, 'distal_point');
        
        
    end
        
    
end