function [nodes edges] = prune_skeleton_to_block(nodes,edges,min_coords, max_coords)
    num_nodes = size(nodes,1);
    num_edges = size(edges,1);

    %first, find interior nodes
    interior_nodes = false(num_nodes,1);
    for n = 1:num_nodes
        interior_nodes(n) = check_inside(nodes(n,:), min_coords, max_coords);
    end
    
    valid_edges = false(num_edges,1);        
    
    
    for n = 1:num_edges
            
        if interior_nodes(edges(n,1)) && interior_nodes(edges(n,2))
            valid_edges(n) = true;
        elseif interior_nodes(edges(n,1)) || interior_nodes(edges(n,2))
            if interior_nodes(edges(n,1))
                good_node = 1;
                bad_node = 2;
            elseif interior_nodes(edges(n,2))
                good_node = 2;
                bad_node = 1;                
            else
            end
            
            gn = nodes(edges(n,good_node),:);
            bn = nodes(edges(n,bad_node),:);
            
            [valid_edges(n) new_node] = rect_search_t(gn, bn, min_coords, max_coords, 12);
            
            if valid_edges(n)
                nodes(end+1,:) = new_node;
                edges(n,bad_node) = size(nodes,1);
            end
        else
            gn = nodes(edges(n,1),:);
            bn = nodes(edges(n,2),:);
            
            [valid_edges(n) new_node] = rect_search_t(gn, bn, min_coords, max_coords, 12);
            if ~valid_edges(n)                
                [valid_edges(n) new_node] = rect_search_t(bn, gn, min_coords, max_coords, 12);
            end
            
            if valid_edges(n)
                [valid_edges(n) newer_node] = rect_search_t(new_node, gn, min_coords, max_coords, 12);
                if ~valid_edges(n)
                    [valid_edges(n) newer_node] = rect_search_t(new_node, bn, min_coords, max_coords, 12);
                end
                
                nodes(end+(1:2),:) = [new_node; newer_node];
                edges(n,:) = size(nodes,1) + (-1:0);                                
            end
                                    
        end
    end
    
    edges(~valid_edges,:) = [];
    
    [nodes edges] = clean_graph(nodes, edges);
end


function is_in = check_inside(n, min_coords, max_coords)
    if all(n >= min_coords) && all(n <= max_coords)
        is_in = true;
    else
        is_in = false;
    end
end

function [t_found new_node] = rect_search_t(gn, bn, min_coords, max_coords, iters)
    t_inc = .5;
    t = .5;

    t_found = false;
    for iter = 1:iters
        t_inc = t_inc/2;
        new_node = gn*t + bn*(1-t);
        if check_inside(new_node, min_coords, max_coords)
            t_found = true;
            t = t-t_inc;
        else
            t = t+t_inc;
        end
    end
end