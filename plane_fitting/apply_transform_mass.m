function p = apply_transform_mass(T,p)

    for n = 1:3
        p(:,n) = p(:,n).*T.res(n);
    end

    p = p*T.global_Q';
    
    num_nodes = size(T.nodes,1);
    
    node_dist = zeros(size(p,1),num_nodes);
    for n = 1:num_nodes
        node_dist(:,n) = sqrt((p(:,2)-T.nodes(n,1)).^2 + ...
        (p(:,3)-T.nodes(n,2)).^2);
    end

    node_dist = T.reach-node_dist;
    node_dist(node_dist<0) = 0;
    
    node_div = 1./sum(node_dist,2);
    z_trans = zeros(1,num_nodes);
    for n = 1:num_nodes
        z_trans(n) = T.translation(n,1);
    end
    
    p(:,1) = p(:,1) + sum(node_dist.*(node_div*z_trans),2);
    
end