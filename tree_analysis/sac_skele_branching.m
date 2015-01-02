function sac_skele_branching

    C = get_constants;
    
    sacs = C.type.sure_off_sac;
    num_sacs = length(sacs);
    
    
    trunk_loc = 16;
    
    for n = 1:num_sacs;
%         cell_dat = cell_data(sacs(n));
%         
%         sac_mid = cell_dat.get_midpoint(true);
        skele = load(['/net/omicfs/home/matthew/stratification/skeletons/s' num2str(sacs(n)) '.mat']);
        skele.edges = full(skele.edges);
        
        
        node_depth = C.f(skele.nodes(:,1));
        
        num_nodes = size(skele.nodes,1);
        
        cell_dat = cell_data(sacs(n));
        soma_pos = cell_dat.get_midpoint;
        dist_from_soma = zeros(num_nodes,1);
        xy_dist_from_soma = zeros(num_nodes,1);
        for d = 1:3
            dist_from_soma = dist_from_soma + (skele.nodes(:,d)-soma_pos(d)).^2;
        end
        for d = 2:3
            xy_dist_from_soma = xy_dist_from_soma + (skele.nodes(:,d)-soma_pos(d)).^2;
        end
        dist_from_soma = sqrt(dist_from_soma);
        xy_dist_from_soma = sqrt(xy_dist_from_soma);
        [~,soma_node] = min(dist_from_soma);
                
        
        %find branches 
        [is_branch is_leaf] = node_summary(skele);
        
        is_critical = is_branch | is_leaf;
        is_critical(node_depth<trunk_loc & xy_dist_from_soma < 10000) = false;
        is_critical(soma_node) = true;
        
        crit_node_list = find(is_critical);
        
        subskele = remove_trivial_skeleton_nodes(skele, crit_node_list, soma_node);
        subsoma_node = find(crit_node_list==soma_node,1,'first');
        subtree = skele2tree(subskele, subsoma_node);
        
        new_bough = subsoma_node;
        bough_nodes = [];
        while length(new_bough)<5            
            bough_nodes = [bough_nodes new_bough];
            newer_boughs = [];
            for k = new_bough
                newer_boughs = [newer_boughs subtree{k}.children];
            end
            new_bough = newer_boughs;
        end
            
        subsubskeles = split_skele(subskele, bough_nodes);
            
        
        
        figure; hold on; c = colormap('Lines');
        plot_skeleton(subskele, 1, [0 0 0]);
        for k = 1:length(subsubskeles)
            plot_skeleton(subsubskeles{k}, 1, c(k,:));
        end
%         figure; plot_skeleton(subskele, 1, [0 0 0]);
%         hold on; scatter(skele.nodes(soma_node,2), skele.nodes(soma_node,3));
%         
%         figure; plot_skeleton(subskele, 3, [0 0 0]);
        
    end
            
            