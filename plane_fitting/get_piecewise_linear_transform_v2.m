function T = get_piecewise_linear_transform_v2(num_rings)

    T.transform = 'original';

    C = get_constants;
    %get global linear transform
%     Q = get_best_rigid_transform;
    
    res = C.res;
    density = .001;

    off_sacs = [70014 70016 70048 70050 70066 70068 70076 70077 70079 ...
        70080 70081 70082 70083 70084 70085 70086 70087 70088 70089 ...
        70090 70093 70094 70095 70096 70097 70098 70099 70100 70101 ...
        70102 70104 70105 70106 70108 70023 70024 70025 70026 70027 ...
        70030 70031 70032 70033 70034 70035];
    
    surface_dir = '../surface_points/';
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
    
    for k = 1:3
        p(:,k) = p(:,k)*res(k);
    end
    
    Q = find_planar_rotation_iterative(p, .001);    
    p = p*Q';
    
    T.global_Q = Q;
    T.center = mean(p);
    
%     p(:,1) = p(:,1)-T.center(1);
%     p(:,2) = p(:,2)-T.center(2);
%     p(:,3) = p(:,3)-T.center(3);
%         
    dist = sqrt((p(:,2)-T.center(1)).^2 + (p(:,3)-T.center(2)).^2);
            
    T.radius = max(dist)/(num_rings+sqrt(3)/2);
    T.reach = T.radius*sqrt(3)/2;
    
    T.num_nodes = 1+sum(6*(1:num_rings));
    
    T.nodes = zeros(T.num_nodes,2);
    T.quaternion = zeros(T.num_nodes,4);
    T.num_points = zeros(T.num_nodes,1);
    T.translation = zeros(T.num_nodes,3);
    T.depth = zeros(T.num_nodes,1);
    
    is_valid = true(T.num_nodes,1);
    
    
    k = 0;
   
    for r = 0:num_rings
        if r==0
            N = 1;
        else
            N = 6*r;
        end
        
        for n = 1:N
            k = k+1;
            
            theta = 2*pi*n/N;
            remainder_angle = mod(theta,pi/3);
            
            if remainder_angle == 0
                h = T.radius*r;
            else
                gamma = pi*2/3 - remainder_angle;
                h = T.radius*r/sin(gamma)*sin(pi/3);            
            end
%             disp(h);
            
            T.nodes(k,:) = h*[cos(theta), sin(theta)] + T.center(2:3);
            
            dist = sqrt((p(:,2)-T.nodes(k,1)).^2 + (p(:,3)-T.nodes(k,2)).^2);
            valid_p = p(dist<=T.reach,:);
            
            T.depth(k) = mean(valid_p(:,1));
            
            if ~isempty(valid_p)
            
            T.num_points(k) = size(valid_p,1);
            
            Q_xaxis = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            
            [Q, ~, used_p] = find_planar_rotation_iterative(valid_p*Q_xaxis, .001);
            
            Q = Q_xaxis*Q*Q_xaxis';
            
            T.quaternion(k,:) = Q2quaternion(Q);
%             [T.theta(k) T.phi(k) T.psi(k)] = Q2angles(Q);
%             
%             
%             Q = angles2Q(T.theta(k), T.phi(k), T.psi(k));
%             valid_p = valid_p*Q';            
%             figure; scatter(valid_p(:,1),valid_p(:,3),1,'b')
            
            used_p = used_p*Q_xaxis';
            used_p = used_p*Q';  
%             figure;
%         scatter(valid_p(:,1),valid_p(:,3), 1); hold on
%         scatter(qp(:,1),qp(:,3), 1, 'r')
            
            
            T.translation(k,:) = -[mean(used_p(:,1)) 0 0];
            
            
%             
%             
%             
%             
%             off_T = T;
%             off_T.nodes((1:end)~=k,:) = Inf;
%             off_T.global_Q = eye(3);
%                 
%             new_p = zeros(size(valid_p)); 
%             for l = 1:size(valid_p,1);
%                 new_p(l,:) = apply_transform(off_T,valid_p(l,:));
%             end
%             figure;
%             scatter(new_p(:,1),new_p(:,3), 1);
%             figure;
%             scatter(valid_p(:,1),valid_p(:,3), 1, 'r')
%             
%             else
%                 is_valid(k) = false; 
%             end
            
        end
           
    end
    
    
    T.nodes(~is_valid,:) = [];
    T.quaternion(~is_valid,:) = [];
    T.num_points(~is_valid,:) = [];
    T.translation(~is_valid,:) = [];
    
%     for k = 1:size(T.nodes,1)
%             
%         
%         dist = sqrt((p(:,2)-T.nodes(k,1)).^2 + (p(:,3)-T.nodes(k,2)).^2);
%         valid_p = p(dist<=T.reach,:);
% 
%         off_T = T;
%         off_T.nodes((1:end)~=k,:) = Inf;
% 
%         new_p = zeros(size(valid_p)); 
%         for l = 1:size(valid_p,1);
%             new_p(l,:) = apply_transform(off_T,valid_p(l,:));
%         end
%         figure;
%         scatter(new_p(:,1),new_p(:,3), 1); hold on
%         scatter(valid_p(:,1),valid_p(:,3), 1, 'r')
%             
%             
% %         end
%     end   
    
    
    
end
    
    
    
    
    
    
    
    