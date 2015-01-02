function density = get_density_grid(cell_num, varargin)
    %compute density of cell by depth, angle, and distance.
    %If use_soma == true, uses the soma as the reference center. Otherwise
    %uses arbor mean
    
    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addOptional('force_recalc', false, @islogical);
    
    p.parse(cell_num, varargin{:});    
    s = p.Results;
    
   
    C = get_constants;
    
    out_fn = [C.strat_dir '/cell_' num2str(s.cell_num) '_' num2str(C.cell_dsmp_fact(1)) '_dsmp.mat'];
    
    if exist(out_fn,'file') && ~s.force_recalc
        load(out_fn)
        return
    end
    
        
    
    fn = [C.point_dir '/cell_' num2str(s.cell_num) '_surface.mat'];
    load(fn);
    
    
    depth = C.f(surface_points(:,1));
    depth = round(depth) - C.strat_x(1) + 1;
    
    
    x = ceil(surface_points(:,2)/C.cell_dsmp_fact(1));
    y = ceil(surface_points(:,3)/C.cell_dsmp_fact(2));

    is_bad = (depth > C.strat_x(end) - C.strat_x(1) + 1) | depth <= 0 | x<=0 | y<=0;
    depth = ceil(depth(~is_bad)/C.cell_dsmp_fact(3));
    x = x(~is_bad);
    y = y(~is_bad);
    
    max_depth = ceil(length(C.strat_x)/C.cell_dsmp_fact(3));
    
    x = x*max_depth + depth;
    
    if isempty(x)
        density = [];
    else
        density = sparse(x,y,ones(length(x),1),max(x),max(y));
    end
        
    save(out_fn, 'density');
end
    
