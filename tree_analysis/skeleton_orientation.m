function orientations = skeleton_orientation(cell_num, num_samples, depth_range, node_distance)

    C = get_constants;
    
    load([C.skele_dir 's' num2str(cell_num), '.mat']);
    
%     figure; scatter(nodes(:,2), nodes(:,3));
    

    d = C.f(nodes(:,1));
    nodes = nodes(d>depth_range(1) & d<depth_range(2),:);
    
    num_nodes = size(nodes,1);
    orientations = zeros(num_samples, 3);
    
        
    for n = 1:num_samples
        r = ceil(rand*num_nodes);
        
        pos = nodes(r,:);
        d = sqrt(sum((nodes - ones(num_nodes,1)*pos).^2,2));
        
        sub_nodes = nodes(d <= node_distance, :);
        
        sub_nodes = sub_nodes - ones(size(sub_nodes,1),1)*pos;
        C = sub_nodes'*sub_nodes;
        [eigv lambda] = eig(C);
        
        orientations(n,:) = eigv(:,1)';
        
    end
    
    num_bins = 36;
    theta = ((1:num_bins)-.5) / num_bins * 2 * pi - pi;
    
    angles = atan2(orientations(:,3), orientations(:,2));    
    angles = angles(abs(orientations(:,1)) > .2, :);
    angles = [angles; angles+pi];
    angles(angles>pi) = angles(angles>pi) - 2*pi;
    h = hist(angles, theta);
    figure; polar([theta theta(1)], [h h(1)]);
    
end
    