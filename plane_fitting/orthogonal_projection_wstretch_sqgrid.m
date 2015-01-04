function T = orthogonal_projection(grid_size)

    T.transform = 'orthostretch';
        
    

    
    C = get_constants;
    %get global linear transform
%     Q = get_best_rigid_transform;
    
    res = C.res;
    T.res = res;
    
    density = .002;

    [grid_X grid_Y] = meshgrid(...
        C.min_xy(2):(C.max_xy(2) - C.min_xy(2))/(grid_size(2)-1):C.max_xy(2), ...
        C.min_xy(1):(C.max_xy(1) - C.min_xy(1))/(grid_size(1)-1):C.max_xy(1));
        
    
    
    off_sacs = C.type.sure_off_sac;
    on_sacs = C.type.on_sac;
    
    p = get_points(off_sacs, density);    
    op = get_points(on_sacs, density);    
    
    for k = 1:3
        p(:,k) = p(:,k)*res(k);
        op(:,k) = op(:,k)*res(k);
    end
    
    Q = find_planar_rotation_iterative(p, .001);    
    p = p*Q';
    op = op*Q';
    
    percentile_cut = .10;
    
    tp = sort(p(:,1));
%     min_p = tp(round(length(tp)*percentile_cut));
    max_p = tp(round(length(tp)*(1-percentile_cut)));
%     p(p(:,1)<min_p | p(:,1)>max_p,:) = [];
    p(p(:,1)>max_p,:) = [];
    
    tp = sort(op(:,1));
    min_op = tp(round(length(tp)*percentile_cut));
%     max_op = tp(round(length(tp)*(1-percentile_cut)));
%     op(op(:,1)<min_op | op(:,1)>max_op,:) = [];
    op(op(:,1)<min_op,:) = [];
    
    
    
    T.global_Q = Q;
    T.center = mean(p);
    
%     p(:,1) = p(:,1)-T.center(1);
%     p(:,2) = p(:,2)-T.center(2);
%     p(:,3) = p(:,3)-T.center(3);
%         
    dist = sqrt((p(:,2)-T.center(1)).^2 + (p(:,3)-T.center(2)).^2);
%     odist = sqrt((op(:,2)-T.center(1)).^2 + (op(:,3)-T.center(2)).^2);
            
%     T.radius = max(dist)/(num_rings+sqrt(3)/2)*1.2; %bit of a hack here
    T.reach = sqrt((grid_X(1,2)-grid_X(1,1))^2 + (grid_Y(2,1)-grid_Y(1,1))^2);
    
    T.num_nodes = numel(grid_X);
    
    T.nodes = zeros(T.num_nodes,2);
    T.num_points = zeros(T.num_nodes,1);
    T.off_chat = zeros(T.num_nodes,1);
    T.on_chat = zeros(T.num_nodes,1);
    T.depth = zeros(T.num_nodes,1);
    
%     is_valid = true(T.num_nodes,1);
    
    
   
    for k = 1:T.num_nodes
            
            
%             theta = 2*pi*n/N;
%             remainder_angle = mod(theta,pi/3);
            
%             if remainder_angle == 0
%                 h = T.radius*r;
%             else
%                 gamma = pi*2/3 - remainder_angle;
%                 h = T.radius*r/sin(gamma)*sin(pi/3);            
%             end
%             disp(h);
            
            T.nodes(k,:) = [grid_Y(k), grid_X(k)];
            
            dist = sqrt((p(:,2)-T.nodes(k,1)).^2 + (p(:,3)-T.nodes(k,2)).^2);
            odist = sqrt((op(:,2)-T.nodes(k,1)).^2 + (op(:,3)-T.nodes(k,2)).^2);
            valid_p = p(dist<=T.reach,:);
            valid_op = op(odist<=T.reach,:);
            
            T.off_chat(k) = median(valid_p(:,1));            
            T.on_chat(k) = median(valid_op(:,1));            
            
            T.num_points(k) = size(valid_p,1);
            
            
            
        
           
    end
    
    figure; scatter(T.nodes(:,1), T.nodes(:,2));
    
    is_valid = ~isnan(T.on_chat) & ~isnan(T.off_chat);
    
    T.nodes(~is_valid,:) = [];
    T.num_points(~is_valid,:) = [];
    T.depth(~is_valid,:) = [];
    T.off_chat(~is_valid,:) = [];
    T.on_chat(~is_valid,:) = [];
    
    T.mean_band_distance = mean(T.off_chat - T.on_chat);
    
    
end
    
function p = get_points(off_sacs, density)

    C = get_constants;

    surface_dir = C.surface_point_dir;
    num_cells = length(off_sacs);
    
    fns = cell(num_cells,1);
    total_points = 0;
    for n = 1:num_cells
        fns{n} = [surface_dir 'cell_' num2str(off_sacs(n)) '_surface.mat'];
        fn_data = whos('-file', fns{n});
        total_points = total_points + ceil(fn_data.size(1)*density);
    end
    
    p = zeros(total_points,3);
    k = 0;
    for n = 1:num_cells
        load(fns{n});
        p_to_add = ceil(size(surface_points,1)*density);
        p(k + (1:p_to_add),:) = surface_points(floor(1:1/density:end),:);
        k = k+p_to_add;
    end
    
end
    
    
    
    