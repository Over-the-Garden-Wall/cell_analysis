function [conns, vericose_nodes] = detect_vericose_contacts(cell_num, max_distance_to_node, vericosity_threshold, min_radius, max_contacts_per_vericosity)
    
    C = get_constants;
    
    if ~exist('max_contacts_per_vericosity', 'var') || isempty(max_contacts_per_vericosity)
        max_contacts_per_vericosity = 0;
    end
    
    try
        load([C.skele_dir 's' num2str(cell_num) '.mat']);
    catch ME
        disp(ME.message)
        conns = [];
        vericose_nodes = [];
        return
    end
    
    
    c_d = cell_data(cell_num);
    all_conns = double(c_d.contacts);
    
    soma_pos = c_d.get_midpoint;
    conn_d_from_soma = sqrt((all_conns(4,:) - soma_pos(2)).^2 + (all_conns(5,:) - soma_pos(3)).^2);
    
    all_conns = all_conns(:,conn_d_from_soma > min_radius);
    
    node_d_from_soma = sqrt((nodes(:,2) - soma_pos(2)).^2 + (nodes(:,3) - soma_pos(3)).^2);
    nodes = nodes(node_d_from_soma > min_radius,:);
    node_diameter = node_diameter(node_d_from_soma > min_radius);

    
    nodes = nodes(node_diameter > vericosity_threshold,:);    
    
    dist_threshold = vericosity_threshold * 50;
        num_conns = size(all_conns,2);
    
    if max_contacts_per_vericosity > 0
        %consolidate nodes, n^2 method
        num_nodes = size(nodes,1);
        node2node_dist = sqrt((nodes(:,1)*ones(1,num_nodes) - ones(num_nodes,1)*nodes(:,1)').^2 + ...
            (nodes(:,2)*ones(1,num_nodes) - ones(num_nodes,1)*nodes(:,2)').^2 + ...
            (nodes(:,3)*ones(1,num_nodes) - ones(num_nodes,1)*nodes(:,3)').^2);


        cluster_nums = zeros(num_nodes,1);
        max_cluster = 0;
        for n = 1:num_nodes
            my_cluster = cluster_nums(n);
            if my_cluster == 0
                max_cluster = max_cluster + 1;
                my_cluster = max_cluster;
            end

            t = node2node_dist(n,:) < dist_threshold;
            cluster_nums(t) = my_cluster;        
        end

        cluster_center = zeros(max_cluster,3);
        for n = 1:max_cluster
            abs_center = mean(nodes(cluster_nums==n,:), 1);
            node_d_from_center = sqrt((nodes(:,1) - abs_center(1)).^2 + ...
                (nodes(:,2) - abs_center(2)).^2 + ...
                (nodes(:,3) - abs_center(3)).^2);
            [dummy, closest_node] = min(node_d_from_center);
            cluster_center(n,:) = nodes(closest_node, :);
        end
        nodes = cluster_center;
    
        num_ver = size(nodes,1);
        d_from_vericosity = sqrt(...
            (nodes(:,2)*ones(1,num_conns) - ones(num_ver,1)*all_conns(4,:)).^2 + ...
            (nodes(:,3)*ones(1,num_conns) - ones(num_ver,1)*all_conns(5,:)).^2);
        
        is_near_vericosity = false(1,num_conns);
        for n = 1:size(nodes,1)
            valid_contacts = find(d_from_vericosity(n,:) <= max_distance_to_node);
            for k = 1:min(max_contacts_per_vericosity, length(valid_contacts));
                [dummy, max_ind] = max(all_conns(2,valid_contacts));
                all_conns(2,max_ind) = 0;
                is_near_vericosity(valid_contacts(max_ind)) = true;
            end
        end
        
    else
                
        num_ver = size(nodes,1);
        d_from_vericosity = sqrt(...
            (nodes(:,2)*ones(1,num_conns) - ones(num_ver,1)*all_conns(4,:)).^2 + ...
            (nodes(:,3)*ones(1,num_conns) - ones(num_ver,1)*all_conns(5,:)).^2);
        is_near_vericosity = any(d_from_vericosity <= max_distance_to_node);
    end

    %     toc
    %     any(i_n_v ~= is_near_vericosity)
    %     detect_vericose_contacts(70014, 4000, 3, 60000);

        conns = all_conns(:, is_near_vericosity);
        vericose_nodes = nodes;
    
end
        