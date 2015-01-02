function create_subsac(cell_num)
    C = get_constants;

    skele = load(['/data/home/greenem/stratification/tree_analysis/skeletons/s' num2str(cell_num) '.mat']);
    skele.edges = full(skele.edges);

    
    fh = figure;
    
    scatter(skele.nodes(:,1),skele.nodes(:,2), 1);
    [z1, x] = getpts(fh);
    
    close(fh);
    
    fh = figure;
    
    scatter(skele.nodes(:,1),skele.nodes(:,3), 1);
    [z2, y] = getpts(fh);
    
    close(fh);
    
    to_remove = skele.nodes(:,2) < x(1) & skele.nodes(:,2) > x(2);
    to_remove = to_remove & skele.nodes(:,3) < y(1) & skele.nodes(:,2) > y(2);
    to_remove = to_remove & skele.nodes(:,3) < mean([z1(:); z2(:)]);
    
    subskeles = split_skele(skele, find(to_remove));
    
    figure; hold all
    
    scatter(skele.nodes(:,2),skele.nodes(:,3), 1);
    for k = 1:length(subskeles)
        scatter(subskeles{k}.nodes(:,2),subskeles{k}.nodes(:,3), 1);
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