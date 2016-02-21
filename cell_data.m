classdef cell_data < handle
    properties
        cell_id = 0;
        V = 0;
        SA = 0;
        has_soma = false;
        is_symmetric = false;
        type = 'unknown';
        subtype = 'unknown';
        
        dist_axis = [0 0];        
        
        percentile_depth = zeros(100,1);
        stratification = [];
        
        contacts = [];
        contact_map = containers.Map;
        contact_area = [];
        contact_count = [];
        contact_ids = [];
        
        hull_2d = [];
        hull_area = 0;
        
    end
    
    
    properties (Access=protected)
        version = -1;
        
        %these are private because we want to force a flag for whether to
        %use mean point or soma point
        meanLoc = [0 0 0];
        somaLoc = [0 0 0];
        
        point_count = [];
        point_count_grid = [];
        
        surface_file = [];
        point_files = {};
        
    end
    
    properties (Access=protected, Hidden)        
        point_count_loaded = false;
        
    end
    
    methods
        
        function cd = cell_data(cell_num, force_construction)
            
            if ~exist('force_construction', 'var') || isempty(force_construction)
                force_construction = false;
            end
            
            
            C = get_constants;
            fn = [C.celldata_dir '/cell_' num2str(cell_num) '_data.mat'];
            if ~exist(fn,'file') || force_construction
                disp(['Constructing cell data for ' num2str(cell_num)]);
            else
                load(fn);
                if cd.version ~= C.celldata_version
                    disp(['Cell data for ' num2str(cell_num) ' is out of date. Updating now']);
                else
                    return
                end
            end
            
            %make data
            
            cd.cell_id = cell_num;
            cd.version = C.celldata_version;
            
            if cell_num < 20000
                cd.type = 'ganglion';
            elseif cell_num < 70000
                cd.type = 'bipolar';
            elseif cell_num < 80000
                cd.type = 'amacrine';
                if cell_num > 71000
                    cd.subtype = 'subarbor';
                else
                    cd.subtype = 'full';
                end
            else
                cd.type = 'ganglion';
            end
            
            fns = fieldnames(C.type);
            for k = 1:length(fns)
                if any(cell_num==C.type.(fns{k}))
                    cd.subtype = fns{k};
                end
            end
            
            [cd.SA cd.V] = get_size_stats(cell_num);
            
            if strcmp(cd.type, 'bipolar');
                cd.has_soma = false;
            else
                cd.has_soma = true;
            end
            
            cd.meanLoc = get_mean_loc(cell_num);
            cd.somaLoc = get_soma_loc(cell_num);
            
            
            
            cd.surface_file = [C.point_dir 'cell_' num2str(cell_num) '_surface.mat'];
            
            if strcmp(cd.type, 'ganglion') || (strcmp(cd.type, 'amacrine') && strcmp(cd.subtype, 'subarbor'));
                cd.is_symmetric = false;
                
                load(cd.surface_file);
                surface_points(C.f(surface_points(1,:)) > 80,:) = [];
                surface_points = double(surface_points(:,2:3));
                for d = 1:2
                    surface_points(:,d) = (surface_points(:,d)-cd.somaLoc(d+1))/size(surface_points,1);
                end

                distal_p = get_distal_loc(cell_num);
                
                daxis = distal_p(2:3) - cd.somaLoc(2:3);
                daxis = daxis/norm(daxis);
                theta = atan2(daxis(2),daxis(1));
                step = pi/8;
                step_dir = 1;
                
                num_p = size(surface_points,1);
                
                while 1
                    orth_d = [-daxis(2) daxis(1)];
                    num_pos = sum(surface_points(:,1)*orth_d(1) + surface_points(:,2)*orth_d(2)>0);
%                     disp(num_pos/num_p)
                    if abs(num_pos/num_p - .5) < .001
                        break
                    elseif num_pos/num_p < .5
                        %too few in pos direction, sweep ccw
                        if step_dir == 1
                            step = step/2;
                            step_dir = -1;
                        end
                    elseif num_pos/num_p > .5
                        if step_dir == -1
                            step = step/2;
                            step_dir = 1;
                        end
                    end
                            
                    theta = theta + step_dir*step;
                    daxis = [cos(theta) sin(theta)];

                end
                cd.dist_axis = daxis;
                
%                 
%                 cov = surface_points'*surface_points;
%                 [vecs vals] = svd(cov);
%                 
%                 cd.dist_axis = vecs(:,1);
                
                                
            else
                cd.is_symmetric = true;
                            
                            
            end
            
            
            pc = get_density_all(cell_num, true);
            cd.stratification = sum(sum(pc,3),2);
            cd.stratification = cd.stratification / sum(cd.stratification);
            
            cum_strat = cumsum(cd.stratification);
            
            delta = .00001;
            for q = 1:100;
                q_val = q/100;
                q_point = find(cum_strat>=q_val-delta,1,'first');
                cd.percentile_depth(q) = C.strat_x(q_point);
            end
            
            c = cd.get_contacts;
            
            cd.contact_ids = unique(c(1,:));
            num_ids = length(cd.contact_ids);
            
            map_data = cell(num_ids,2);
            
            for k = 1:num_ids;
                map_data{k,1} = k;
                map_data{k,2} = cd.contact_ids(k);
            end
            
            if isempty(map_data)
                cd.contact_map = [];
            else
                cd.contact_map = containers.Map(map_data(:,2), map_data(:,1));
            end
            cd.contact_area = zeros(1,num_ids);
            cd.contact_count = zeros(1,num_ids);
                        
            for k = 1:num_ids;
                is_me = c(1,:) == cd.contact_ids(k);
                cd.contact_area(k) = sum(c(2,is_me));
                cd.contact_count(k) = sum(is_me);
            end
            
            cd.contacts = c;
                        
            cd.point_files = get_files_with_names_including(C.raw_point_dir,num2str(cell_num));
            for k = 1:length(cd.point_files)
                cd.point_files{k} = [C.raw_point_dir '/' cd.point_files{k}];
            end
            
%             if strcmp(cd.type,'bipolar')
                cd.hull_2d = cd.get_2d_hull(1);
                
                cd.hull_area = polyarea(cd.hull_2d(:,1), cd.hull_2d(:,2));
                
%             else
%                 cd.hull_2d = [];
%                 cd.hull_area = [];
%             end
            
            save(fn,'cd');
        end
        
        function value = get_point_count(cd)
            if ~cd.point_count_loaded
                cd.point_count = get_density_all(cd.cell_id);           
                cd.point_count_loaded = true;
            end
            value = cd.point_count;
        end
        
        function value = get_point_count_grid(cd)
            if isempty(cd.point_count_grid)
                cd.point_count_grid = get_density_grid(cd.cell_id);           
            end
            value = cd.point_count_grid;
        end
          
    
        function midpoint = get_midpoint(cd, use_soma)
            if ~exist('use_soma','var') || isempty(use_soma)
                use_soma = cd.has_soma;
            end

            if use_soma
                midpoint = cd.somaLoc;
            else
                midpoint = cd.meanLoc;
            end
        end
        
        function c = get_contacts(cd)
            if isempty(cd.contacts)
                C = get_constants;
                load(C.conn_loc);
                
                latter_is_target = conns(2,:) == cd.cell_id;
                former_is_target = conns(1,:) == cd.cell_id;
    
                conns(1:2,latter_is_target) = conns([2 1], latter_is_target);
                conns = conns(:, latter_is_target | former_is_target);
    
                cd.contacts = conns(2:6,:);
            end
            c = cd.contacts;           
        end
        
        function surface_points = get_surface(cd, bounds)            
            load(cd.surface_file);
            
            if exist('bounds','var') && ~isempty(bounds)
                for d = 1:3
                    is_valid = ...
                        surface_points(:,d) >= bounds(d,1) || ...
                        surface_points(:,d) <= bounds(d,2);
                    surface_points = surface_points(is_valid,:);
                end                
            end
            
        end
        
        
        function vol_points = get_volume(cd, bounds)
            if ~exist('bounds','var') || ~all(size(bounds)==[3 2])
                error('get_volume requires bounds argument of size [3 2]')
            end
            bounds = uint32(bounds);
            disp('getting points! (might be slow)'); tic
            
            vol_points = [];
            for k = 1:length(cd.point_files);
                load(cd.point_files{k});
                for d = 1:3
                    is_valid = ...
                        p(:,d) >= bounds(d,1) & ...
                        p(:,d) <= bounds(d,2);
                    p = p(is_valid,:);
                end     
                vol_points = [vol_points; p];
            end
            
            toc
            
        end
        
        function bounds = get_volume_bounds(cd)            
            disp('getting points! (might be slow)'); tic
            
            bounds = [Inf 0; Inf 0; Inf 0];
            
            for k = 1:length(cd.point_files);
                load(cd.point_files{k});
                bounds(:,1) = min(bounds(:,1), double(min(p)'));
                bounds(:,2) = max(bounds(:,2), double(max(p)'));                
            end
            
            toc
            
        end
        
        function hull = get_2d_hull(cd, is_on_construction)
            if ~exist('is_on_construction','var')
                warning('get_2d_hull is deprecated, use hull_2d property');
            end
            
            p = cd.get_surface;
%             
%             region_size = 3000;
%             
%             p = p(:,2:3);
%             minp = min(p);
%             maxp = max(p);
%             
%             [X Y] = meshgrid(minp(1)-region_size/4:region_size/2:maxp(1)+region_size/4, ... 
%                 minp(2)-region_size/4:region_size/2:maxp(2)+region_size/4);
%             
%             hull = [];
%             for x = 1:size(X,1)
%                 for y = 1:size(X,2);
%                     sub_p = p(abs(p(:,1)-X(x,y)) <= region_size/2 & abs(p(:,2)-Y(x,y)) <= region_size/2,:);
%                     if ~isempty(sub_p) && size(sub_p,1) > 2
%                         try
%                             k = convhull(sub_p(:,1),sub_p(:,2));
%                             subhull = [];
%                             [subhull(:,1) subhull(:,2)] = poly2cw(sub_p(k,1), sub_p(k,2));
%                             if isempty(hull)
%                                 hull = subhull;
%                             else
%                                 temphull = [];
%                                 [temphull(:,1) temphull(:,2)] = polybool('union', subhull(:,1), subhull(:,2), hull(:,1), hull(:,2));
%                                 hull = temphull;
%                             end
%                         catch ME
%                             disp(ME.message);
%                         end
%                     end
%                 end
%             end
%             
%             % clean internal bits and disconnected regions (which aren't
%             % supposed to exist
%             
%             nan_locs = [0; find(isnan(hull(:,1))); size(hull,1)+1];
%             num_pieces = length(nan_locs)-1;
%             
%             is_inside = false(num_pieces);
%             piece_size = zeros(num_pieces,1);
%             
%             for k = 1:num_pieces
%                 rk = nan_locs(k)+1:nan_locs(k+1)-1;
%                 for j = 1:num_pieces
%                     if k ~= j
%                         rj = nan_locs(j)+1:nan_locs(j+1)-1;                
%                         is_inside(k,j) = inpolygon(hull(rj(1),1), hull(rj(1),2), hull(rk,1), hull(rk,2));
%                     end
%                 end
%                 
%                 piece_size(k) = poly_area(hull(rk,:));
%             end
%             
%             external_pieces = find(~any(is_inside));
%             
%             [max_size max_ind] = max(piece_size(external_pieces));
%             pick = external_pieces(max_ind);
%             
%             hull = hull(nan_locs(pick)+1:nan_locs(pick+1)-1,:);
%             
%             
            
            k = convhull(p(:,2),p(:,3));
            hull = [p(k,2), p(k,3)];
        end
        
        function hull = get_3d_hull(cd)
            p = cd.get_surface;
            k = convhulln(p);
            hull = p(k,:);
        end
        
        
        
    end
end

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

function [SA, vol] = get_size_stats(cell_nums)

    C = get_constants;

    vol_dir = C.raw_point_dir;
    
    sa_dir = C.point_dir;

    num_cells = length(cell_nums);
    SA = zeros(num_cells,1);
    vol = zeros(num_cells,1);
    
    for n = 1:num_cells

        fn =[sa_dir '/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            point_info = whos('-file', fn);
            SA(n) = point_info.size(1);

            fns = get_files_with_names_including(vol_dir, num2str(cell_nums(n)));

            for k = 1:length(fns)
                point_info = whos('-file', [vol_dir fns{k}]);
                vol(n) = vol(n)+point_info.size(1);
            end



        end
    
    end
end

    
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

function mean_point = get_mean_point(cell_num, use_soma)

    if use_soma
        mean_point = get_soma_loc(cell_num); 
    else
        mean_point = get_mean_loc(cell_num);
    end
    
end
