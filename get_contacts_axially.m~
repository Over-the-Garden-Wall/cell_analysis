function [density density_counts] = get_contacts_axially(cell_num, contact_cells, varargin)
    %compute density of cell by depth, angle, and distance.
    %If use_soma == true, uses the soma as the reference center. Otherwise
    %uses arbor mean
    
    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addRequired('contact_cells', @isnumeric);
    p.addOptional('use_soma', check_to_use_soma(cell_num), @islogical);
    
    p.parse(cell_num, contact_cells, varargin{:});    
    s = p.Results;
    
    C = get_constants;
    
    mean_point = get_mean_point(s.cell_num, s.use_soma); 
    
    load(C.conn_loc);
    
    sub_conns = double(pick_conns(conns, s.cell_num, s.contact_cells, false));
    
    
    
        num_points = size(sub_conns,2);


        distal_point = get_distal_loc(s.cell_num);
        
        sub_conns(5,:) = sub_conns(5,:)-mean_point(2);
        sub_conns(6,:) = sub_conns(6,:)-mean_point(3);
        
        dist_vec = distal_point(2:3) - mean_point(2:3);
        dist_vec = dist_vec / sqrt(sum(dist_vec.^2));
        
        dist = sub_conns(5,:)*dist_vec(1) + sub_conns(6,:)*dist_vec(2);
        depth = C.f(sub_conns(4,:));
        angle = atan2(sub_conns(5,:) - mean_point(2), sub_conns(6,:) - mean_point(3));

        depth = round(depth) - C.strat_x(1) + 1;
        dist = ceil(dist/1000) - C.axial_x_min + 1;
        angle = ceil((angle+pi)/C.angle_step);


        is_bad = (depth > C.strat_x(end) - C.strat_x(1) + 1) | depth <= 0 | dist <= 0;
        

        depth(is_bad) = [];
        dist(is_bad) =  [];
        angle(is_bad) = [];
        
        density = zeros(C.strat_x(end) - C.strat_x(1) + 1, max(dist), ceil(2*pi/C.angle_step));
        density_counts = zeros(size(density));
        
        num_points = length(depth);

        for k = 1:num_points
            density(depth(k), dist(k), angle(k)) = density(depth(k), dist(k), angle(k))+sub_conns(3,k);
            density_counts(depth(k), dist(k), angle(k)) = density_counts(depth(k), dist(k), angle(k))+1;            
        end    
    
    
end
    
