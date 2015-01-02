function s = sac_contacts_subsections(cell_nums, cont_types)
    C = get_constants;

    
    s.all_branch_loc = [];
    s.all_is_leaf = [];
    s.all_num_contacts = [];
    s.all_total_contact = [];
    s.branch_num_points = [];
    
    
    num_types = length(cont_types);
    
    hull_intersection = [];
    for k = 1:num_types;
        type_hull = [];
        for c = C.type.(cont_types{k});
            c_d = cell_data(c);
            h = [];
            [h(:,1), h(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
            temp_hull = [];
            if isempty(type_hull)
                type_hull = h;
            else
                [temp_hull(:,1), temp_hull(:,2)] = polybool('union', type_hull(:,1), type_hull(:,2), h(:,1), h(:,2));
                type_hull = temp_hull;
            end
        end
        if k == 1
            hull_intersection = type_hull;
        else
            temp_hull = [];
            [temp_hull(:,1), temp_hull(:,2)] = polybool('intersection', type_hull(:,1), type_hull(:,2), hull_intersection(:,1), hull_intersection(:,2));
            hull_intersection = temp_hull;
        end
            
    end
    
    
    for cell_num = cell_nums
    
        skele = load(['/data/home/greenem/stratification/tree_analysis/skeletons/s' num2str(cell_num) '.mat']);
        skele.edges = full(skele.edges);

        G = skeleton2graph(skele);
        
        num_nodes = size(skele.nodes,1);

        cell_dat = cell_data(cell_num);
        soma_pos = cell_dat.get_midpoint(true);
        dist_from_soma = zeros(num_nodes,1);
        for d = 1:3
            dist_from_soma = dist_from_soma + (skele.nodes(:,d)-soma_pos(d)).^2;
        end
        dist_from_soma = sqrt(dist_from_soma);
        [~,soma_node] = min(dist_from_soma);

        %find branches 
        [is_branch is_leaf] = node_summary(skele);

        
        is_critical = is_branch | is_leaf;
        is_critical(soma_node) = true;

        crit_node_list = find(is_critical);

        subskele = remove_trivial_skeleton_nodes(skele, crit_node_list, soma_node);
        num_subnodes = size(subskele.nodes,1);
        
        is_leaf = false(num_subnodes,1);
        for k = 1:num_subnodes
            is_leaf(k) = sum(subskele.edges(:)==k) == 1;
        end
        
        num_branches = size(subskele.edges,1);
        is_leaf_edge = is_leaf(subskele.edges(:,1)) | is_leaf(subskele.edges(:,2));
        
        
        conts = double(cell_dat.contacts);
        is_valid = inpolygon(conts(4,:)', conts(5,:)', hull_intersection(:,1), hull_intersection(:,2));
        conts = conts(:,is_valid);
        
        num_conts = size(conts,2);
        
        
        cont_type = zeros(1, num_conts);
        
        
        for k = 1:num_types
            for n = 1:num_conts
            
                if any(conts(1,n)==C.type.(cont_types{k}))
                    cont_type(n) = k;
                end
            end
        end
        
        conts(:,cont_type == 0) = [];
        num_conts = size(conts,2);
        cont_type(cont_type == 0) = [];
        
        nearest_node = zeros(num_conts,1);
        for k = 1:num_conts
            d = zeros(num_nodes,1);
            for dim = 1:3
                d = d + (skele.nodes(:,dim) - conts(dim+2,k)).^2;
            end
            [dummy, nearest_node(k)] = min(d);
        end
       
        
        membership = zeros(num_nodes,1);
        
        branch_loc = zeros(num_branches,3);
        branch_num_contacts = zeros(num_branches,num_types);
        branch_total_contact = zeros(num_branches,num_types);
        branch_num_points = zeros(num_branches,1);
        
        
%         figure; hold on
        for n = 1:num_branches;
            [dummy, spath] = graphshortestpath(G,crit_node_list(subskele.edges(n,1)), crit_node_list(subskele.edges(n,2)));
            
%             scatter(skele.nodes(spath,2), skele.nodes(spath,3), 1);
            
            membership(spath) = n;
            branch_loc(n,:) = mean(skele.nodes(spath,:));
            for d = 2:3
                branch_loc(n,d) = branch_loc(n,d) - soma_pos(d);
            end
            
            branch_num_points(n) = length(spath);
        end
        
        for n = 1:num_conts;
            my_type = cont_type(n);
            my_branch = membership(nearest_node(n));
            
            branch_num_contacts(my_branch, my_type) = branch_num_contacts(my_branch, my_type) + 1;
            branch_total_contact(my_branch, my_type) = branch_total_contact(my_branch, my_type) + conts(2,n);
            
        end
        
        max_branch_total_density = max(sum(branch_total_contact,2)./branch_num_points);
        
        figure; hold on
        for n = 1:num_branches;
            h = branch_total_contact(n,1) / sum(branch_total_contact(n,:));
            if isnan(h)
                c = [.8 .8 .8];
            else
%                 w = sum(branch_total_contact(n,:)/branch_num_points(n))/max_branch_total_density;
                w = branch_num_contacts(n) > 3;
                c = (h*[1 1 0] + (1-h)*[0 0 1]) * w + (1-w)*[.8 .8 .8];
            end
            
            [dummy, spath] = graphshortestpath(G,crit_node_list(subskele.edges(n,1)), crit_node_list(subskele.edges(n,2)));
                        
            scatter(skele.nodes(spath,2), skele.nodes(spath,3), 10, 'markerEdgeColor', c);
        end
        
        
        
        s.all_branch_loc = [s.all_branch_loc; branch_loc];
        s.all_is_leaf = [s.all_is_leaf; is_leaf_edge];
        s.all_num_contacts = [s.all_num_contacts; branch_num_contacts];
        s.all_total_contact = [s.all_total_contact; branch_total_contact];
        s.branch_num_points = [s.branch_num_points; branch_num_points];
        
    end
    
end


function near_nodes = nearest_nodes_to_clicks(fh, subskele)
    [x, y] = getpts(fh);
    near_nodes = zeros(length(x),1);
    for n = 1:length(x)
        d = sqrt((subskele.nodes(:,2)-x(n)).^2 + (subskele.nodes(:,3)-y(n)).^2);
        [dummy, near_nodes(n)] = min(d);
    end
end
    
function skele = combine_skeletons(skels)
    %stupid way for now
    num_skels = length(skels);
    
    skele.nodes = zeros(0,3);
    skele.edges = zeros(0,2);
    for n = 1:num_skels
        k = size(skele.nodes,1);
        skele.nodes = [skele.nodes; skels{n}.nodes];        
        skele.edges = [skele.edges; skels{n}.edges + k];
    end
    
end

function [new_skele all_p] = minimal_skele_including(skele, needed_nodes)

    num_needed = size(needed_nodes,1);
    num_nodes = size(skele.nodes,1);
    is_needed = false(num_nodes,1);
    
    for n = 1:num_needed
        k = find(skele.nodes(:,1) == needed_nodes(n,1) & ...
            skele.nodes(:,2) == needed_nodes(n,2) & ...
            skele.nodes(:,3) == needed_nodes(n,3),1,'first');
        is_needed(k) = true;
    end

    G = skeleton2graph(skele);
    [d, p] = graphshortestpath(G, find(is_needed,1,'first'));
    p = p(is_needed);
    all_p = [];
    for n = 1:length(p)
        all_p = [all_p, p{n}];
    end
    all_p = unique(all_p);
    
    nodes2new_nodes = zeros(num_nodes,1);
    for n = 1:length(all_p)
        nodes2new_nodes(all_p(n)) = n;
    end
    
    new_skele.nodes = skele.nodes(all_p,:);
    new_skele.edges = nodes2new_nodes(skele.edges);
    new_skele.edges(new_skele.edges(:,1)==0 | new_skele.edges(:,2)==0,:) = [];
    
    
end