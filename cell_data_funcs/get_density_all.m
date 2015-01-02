function density = get_density_all(cell_num, varargin)
    %compute density of cell by depth, angle, and distance.
    %If use_soma == true, uses the soma as the reference center. Otherwise
    %uses arbor mean
    
    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addOptional('force_recalc', false, @islogical);
    p.addOptional('use_soma', check_to_use_soma(cell_num), @islogical);
    
    p.parse(cell_num, varargin{:});    
    s = p.Results;
    
   
    C = get_constants;
    
    if s.use_soma
        out_fn = [C.strat_dir '/cell_' num2str(s.cell_num) '_soma_' 'all.mat'];
    else
        out_fn = [C.strat_dir '/cell_' num2str(s.cell_num) '_arbor_' 'all.mat'];
    end
    
    if exist(out_fn,'file') && ~s.force_recalc
        load(out_fn)
        return
    end
    
    
    
    
        
    
    mean_point = get_mean_point(s.cell_num, s.use_soma); 
    

    fn = [C.point_dir '/cell_' num2str(s.cell_num) '_surface.mat'];
    load(fn);
    
    num_points = size(surface_points,1);
    
    
    dist = sqrt(sum((surface_points(:,2:3) - ones(num_points,1)*mean_point(2:3)).^2,2));
    depth = C.f(surface_points(:,1));
    angle = atan2(surface_points(:,2) - mean_point(2), surface_points(:,3) - mean_point(3));
    
    
    depth = round(depth) - C.strat_x(1) + 1;
    dist = ceil(dist/1000);
    angle = ceil((angle+pi)/C.angle_step);
    
    density = zeros(max(depth), max(dist), max(angle));
    
    
    is_bad = (depth > C.strat_x(end) - C.strat_x(1) + 1) | depth <= 0;
    
    depth(is_bad) = [];
    dist(is_bad) =  [];
    angle(is_bad) = [];
    
    num_points = length(depth);
    
    for k = 1:num_points
        density(depth(k), dist(k), angle(k)) = density(depth(k), dist(k), angle(k))+1;
    end
        
    save(out_fn, 'density');
end
    
