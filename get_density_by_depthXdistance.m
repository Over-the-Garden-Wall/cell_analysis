function [density is_valid] = get_density_by_depth(cell_num, use_soma, dist_mode, force_recalc)


    if ~exist('force_recalc','var') || isempty(force_recalc)
        force_recalc = false;
    end
    
    if ~exist('use_soma','var') || isempty(use_soma)
        use_soma = false;
    end  
    
    
    if ~exist('dist_mode','var') || (~strcmp(dist_mode, 'axis') && ~strcmp(dist_mode, 'point') )
        warning('unknown dist_mode, defaulting to ''point''');
        dist_mode = 'point';
    end
    
    if strcmp(dist_mode, 'axis')        
        dist_coords = [2 3];
    else
        dist_coords = [1 2 3];
    end
    
    C = get_constants;
    
    annuli_fn = [C.strat_dir '/cell_' num2str(cell_num) '_annuli.mat'];
    if exist(annuli_fn,'file') && ~force_recalc
        load(annuli_fn);
        is_valid = true;
    else    
    
        fn = [C.point_dir '/cell_' num2str(cell_num) '_annuli.mat'];
        if exist(fn,'file')
            data = load(fn);
            fields = fieldnames(data);
            p = f(data.(fields{1})(:,1));            
            is_valid = true;
        else
            warning(['data not found for ' fn]);
            density = [];        
            is_valid = false;
            return
        end
        
        if use_soma
            mid_point = get_soma_loc(cell_num);
        else
            mid_point = mean(p);
        end
        num_points = size(p,1);
        
        dist = sqrt(sum((p(:,dist_coords) - ones(num_points,1)*mid_point(dist_coords)).^2, 2));
        dist = ceil(dist/1000);
        dist(dist==0) = [];
        density = zeros(1,max_dist);
        
        for k = 1:num_points
            density(dist(k)) = density(dist(k))+1;
        end
        
        save(annuli_fn, 'density', 'mid_point');
    end
    
    

    if ~exist('force_recalc','var') || isempty(force_recalc)
        force_recalc = false;
    end

    C = get_constants;
    
    strat_fn = [C.strat_dir '/cell_' num2str(cell_num) '_strat.mat'];
    if exist(strat_fn,'file') && ~force_recalc
        load(strat_fn);
        is_valid = true;
    else    
    
        x = C.strat_x;
        
        fn = [C.point_dir '/cell_' num2str(cell_num) '_surface.mat'];
        if exist(fn,'file')
            data = load(fn);
            fields = fieldnames(data);
            p = C.f(data.(fields{1})(:,1));
            density = hist(p,x);            
            is_valid = true;
            save(strat_fn, 'density');
        else
            warning(['data not found for ' fn]);
            density = zeros(size(x));        
            is_valid(n) = false;
        end

    end
    
end