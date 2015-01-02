function T = get_piecewise_linear_transform(num_rings)

    error('this function is deprecated, use get_piecewise_linear_transform_v2');
    %get global linear transform
%     Q = get_best_rigid_transform;
    

    off_sacs = [70014 70016 70048 70050 70066 70068 70076 70077 70079 ...
        70080 70081 70082 70083 70084 70085 70086 70087 70088 70089 ...
        70090 70093 70094 70095 70096 70097 70098 70099 70100 70101 ...
        70102 70104 70105 70106 70108 70023 70024 70025 70026 70027 ...
        70030 70031 70032 70033 70034 70035];
    
    

    on_cells = [4,5,6,10,11,12,15,16,17,27,33];

    dirs = dir('./');
    
%     k = 0;
%     for n = 1:length(dirs)
%         if length(dirs(n).name)>2 && strcmp(dirs(n).name(1:3),'sac')
%             k=k+1;
%         end
%     end
%     
%     data = cell(k,1);
%     k = 0;
%     for n = 1:length(dirs)
%         if length(dirs(n).name)>2 && strcmp(dirs(n).name(1:3),'sac')
%             k=k+1;
%             [data{k}.nodes data{k}.edges] = get_skeleton(dirs(n).name);
%             data{k}.length = get_skele_length(data{k}.nodes, data{k}.edges);
%             data{k}.sac_num = str2double(dirs(n).name(4:end));
%             if any(on_cells==data{k}.sac_num)
%                 data{k}.on_cell = true;
%             else
%                 data{k}.on_cell = false;
%             end
%             
% %             data{k}.nodes = data{k}.nodes*Q';
%         end
%     end
%     
%     
%     %collect point clouds for off cells
%     point_density = 1000;
%     
%     num_points = 0;
%     for n = 1:length(data);
%         if ~data{n}.on_cell
%             num_points = num_points + floor(data{n}.length/point_density);
%         end
%     end
%     
%     k = 0;
%     p = zeros(num_points,3);
%     for n = 1:length(data);
%         if ~data{n}.on_cell
%             my_num_points = floor(data{n}.length/point_density);
%             p(k+(1:my_num_points),:) = tree2points(data{n}.nodes,data{n}.edges,point_density);
%             k = k+my_num_points;
%         end
%     end
%     
%     %
% %     Q = eye(3);
% % p = randn(100000,3); p(:,1) = p(:,1)/100;
% % p(:,1) = .7 + p(:,1) + ((p(:,2)-.5).^2+(p(:,3)-.2).^2)/10;
% % p(p(:,2).^2+p(:,3).^2 > 2,:) = [];
%     %
%     
%     
%     Q = find_planar_rotation_iterative(p, .00001);    
    
    
    
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
            
            T.nodes(k,:) = h*[cos(theta), sin(theta)];
            
            dist = sqrt((p(:,2)-T.nodes(k,1)).^2 + (p(:,3)-T.nodes(k,2)).^2);
            valid_p = p(dist<=T.reach,:);
            T.num_points(k) = size(valid_p,1);
            
            Q_xaxis = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
            
            Q = find_planar_rotation_iterative(valid_p*Q_xaxis, .001);
            
            Q = Q_xaxis*Q*Q_xaxis';
            
            T.quaternion(k,:) = Q2quaternion(Q);
%             [T.theta(k) T.phi(k) T.psi(k)] = Q2angles(Q);
%             
%             
%             Q = angles2Q(T.theta(k), T.phi(k), T.psi(k));
%             valid_p = valid_p*Q';            
%             figure; scatter(valid_p(:,1),valid_p(:,3),1,'b')
            
            
            rot_p = valid_p*Q';  
%             figure;
%         scatter(valid_p(:,1),valid_p(:,3), 1); hold on
%         scatter(qp(:,1),qp(:,3), 1, 'r')
            
            
            T.translation(k,:) = -[mean(rot_p(:,1)) 0 0];
            
            
            
            
            
            
            off_T = T;
            off_T.nodes((1:end)~=k,:) = Inf;
                
            new_p = zeros(size(valid_p)); 
            for l = 1:size(valid_p,1);
                new_p(l,:) = apply_transform(off_T,valid_p(l,:));
            end
            figure;
            scatter(new_p(:,1),new_p(:,3), 1);
            figure;
            scatter(valid_p(:,1),valid_p(:,3), 1, 'r')
            
            
            
        end
           
    end
    
    
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
    
    
    
    
    
    
    
    