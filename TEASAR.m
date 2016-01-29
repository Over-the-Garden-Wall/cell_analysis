function [nodes, edges, root_node, node_radii] = TEASAR(im_points, dsmp_resolution, scale, const)
    %im_points - n x 3 array of point coordinates
    %input_resolution - self explanatory, difference in real length between
    %                   point [1 1 1] and [0 0 0]
    %dsmp_resolution - resolution of the coordinates after downsampling
    %scale and const - "scale" and "const" parameter from TEASAR paper,
    %                  larger values mean fewer hairs
    %
       
    NDIMS = 3;
    MAX_CHUNK_SIZE = 251 * [1 1 1];
    OVERLAP = 9 * [1 1 1];
    
    
    if ~exist('dsmp_resolution','var') || isempty(dsmp_resolution)
        dsmp_resolution = ones(1, NDIMS);
    end    
    if ~exist('scale','var') || isempty(scale)
        scale = 6 / max(dsmp_resolution);
    end
    if ~exist('const','var') || isempty(const)
        const = 6 / max(dsmp_resolution);
    end
    
    
    
    
    for d = 1:NDIMS
        im_points(:,d) = ceil(im_points(:,d)./dsmp_resolution(d));
    end
    
    im_origin = im_points(1,:);
    
    im_max = double(max(im_points))+1;
    max_feasible_cubes = ceil((im_max+1)./MAX_CHUNK_SIZE);
    [cube_coord(:,1), cube_coord(:,2), cube_coord(:,3)] = ind2sub( ...
        max_feasible_cubes*2+1, 1:prod(max_feasible_cubes*2+1));
    for d = 1:3
        cube_coord(:,d) = cube_coord(:,d) - max_feasible_cubes(d) - 1;
    end
    
    find_cube = @(c,v) find(all(c == ones(size(c,1),1)*v,2));
    
    processing_order = zeros(size(cube_coord,1),1);
    processing_order(find_cube(cube_coord, [0 0 0])) = 1;
    
    
    
    
    cube_to_process = 1;
    processing_max = 1;
    while 1
        loc_to_process = cube_coord(find(processing_order==cube_to_process,1,'first'),:);
        
        if isempty(loc_to_process)% || cube_to_process == 3;
            
            break
        end
        disp([cube_to_process loc_to_process])
        
        
        valid_p = true(size(im_points,1),1);
        for d = 1:3            
            im_bounds(d,1) = im_origin(d) + loc_to_process(d) * (MAX_CHUNK_SIZE(d) - OVERLAP(d)) - MAX_CHUNK_SIZE(d)/2;
            im_bounds(d,2) = im_bounds(d,1) + MAX_CHUNK_SIZE(d) - 1;
            
            valid_p = valid_p & (im_points(:,d) >= im_bounds(d,1) & im_points(:,d) <= im_bounds(d,2));
        end
        my_points = double(im_points(valid_p,:));
        
        for d = 1:3
            my_points(:,d) = my_points(:,d) - double(im_bounds(d,1)) + 1;
        end
        
        d = double(im_bounds);

        
        tic
        cube_out{cube_to_process} = cube_TEASAR(my_points, scale, const);
        cube_out{cube_to_process}.bounds = im_bounds;
        cube_out{cube_to_process}.location = loc_to_process;
        toc
        
%         figure; plot_skeleton(cube_out{cube_to_process},1,[0 0 0]);
        
        %check the six faces for new cubes
        cube_neighborhood = [-eye(3); eye(3)];
        cube_neighborhood_reldim = [(1:3)' ones(3,2); (1:3)' ones(3,1)*2 MAX_CHUNK_SIZE'];
        point_summary = [min(my_points); max(my_points)];
        
        for d = 1:size(cube_neighborhood,1)
            if point_summary(cube_neighborhood_reldim(d,2),cube_neighborhood_reldim(d,1))==cube_neighborhood_reldim(d,3) 
                t = processing_order(find_cube(cube_coord, loc_to_process + cube_neighborhood(d,:)));
                if ~isempty(t) && t == 0
                    processing_max = processing_max+1;
                    processing_order(find_cube(cube_coord, loc_to_process + cube_neighborhood(d,:))) = processing_max;
                end
            end
        end
        cube_to_process = cube_to_process+1;
    end
    %stitch cubes
    
    num_cubes = length(cube_out);
    num_cube_points = zeros(num_cubes,1);
    num_cube_edges = zeros(num_cubes,1);
    for src_cube = 1:num_cubes       
        num_cube_points(src_cube) = length(cube_out{src_cube}.node_radii);
        num_cube_edges(src_cube) = size(cube_out{src_cube}.edges,1);
    end
    num_cube_points_end = cumsum(num_cube_points);
    num_cube_points_start = 1 + [0; num_cube_points_end(1:end-1)];
    num_cube_edges_end = cumsum(num_cube_edges);
    num_cube_edges_start = 1 + [0; num_cube_edges_end(1:end-1)];
    
    num_nodes = num_cube_points_end(end);
    num_edges = num_cube_edges_end(end);
    
    
    nodes = zeros(num_nodes,3);
    node_radii = zeros(num_nodes,1);
    root_node = 1;
    edges = zeros(num_edges,2);
        
    for src_cube = 1:num_cubes
        n = double(cube_out{src_cube}.nodes);
        for d = 1:3
            n(:,d) = n(:,d) + double(cube_out{src_cube}.bounds(d,1)) - 1;
            n(:,d) = n(:,d) * dsmp_resolution(d);
        end
        nodes(num_cube_points_start(src_cube):num_cube_points_end(src_cube), :) = n;
        node_radii(num_cube_points_start(src_cube):num_cube_points_end(src_cube)) = cube_out{src_cube}.node_radii;        
        
        edges(num_cube_edges_start(src_cube):num_cube_edges_end(src_cube), :) = cube_out{src_cube}.edges - 1 + num_cube_points_start(src_cube);
    end
    
    disp('resolving edges'); tic;
    dup_nodes = find_duplicate_entries(nodes);    
    for n = size(dup_nodes,1):-1:1;
        edges(edges==dup_nodes(n,2)) = dup_nodes(n,1);
    end
    node_list = unique(edges(:));
    
    nodes = nodes(node_list,:);
    node_radii = node_radii(node_list);
    for n = 1:length(node_list)
        edges(edges==node_list(n)) = n;
    end
    edges = unique(edges, 'rows');
    toc
    
    %this is a hack to fix an uncommon issue that I haven't identified
    edge_len = sqrt(sum((nodes(edges(:,1),:) - nodes(edges(:,2),:)).^2, 2));
    edges(edge_len > 10,:) = [];
        
end

function cube_out = cube_TEASAR(p, scale, const)
    NDIMS = size(p,2);


    quick_ind = @(S, I) (sub2ind(S, I(:,1), I(:,2), I(:,3)));
    im_max = max(p);
    
    bin_im = false(im_max);
    
    bin_im(quick_ind(im_max, p)) = true;
    
    node2ind = find(bin_im(:));
    
    loc2node = zeros(size(bin_im));
    loc2node(node2ind) = 1:length(node2ind);
    
    N = length(node2ind);
    
    
    ind2node = sparse(node2ind,ones(N,1),(1:N)',prod(im_max),1);
    
    node2loc = zeros(N, 3);
    [node2loc(:,1), node2loc(:,2), node2loc(:,3)] = ind2sub(im_max, node2ind); 
    
    
    
    DBF = bwdist(~bin_im);
    
    
    %below is a change I added to avoid problems with flat portions of the
    %neurite. Change DBF_adjustment_constant to 0 or delete to remove.
    
    DBF_adjustment_constant = (sqrt(2)-1)/100; %pretty small number, this is just to break ties
    for d = 1:NDIMS
        not_d = 1:NDIMS;
        not_d(d) = [];
        
        temp_bin_im = permute(bin_im, [not_d d]);
        DBF_inc = zeros(size(temp_bin_im));
        for n = 1:size(bin_im,d)
            DBF_inc(:,:,n) = bwdist(~temp_bin_im(:,:,n));
        end
        
        DBF = DBF + DBF_adjustment_constant*permute(DBF_inc, [1:d-1 NDIMS d:NDIMS-1]);
    end
            
            
    
    
    M = max(DBF(:))*1.01;
    
    p_v = 5000*(1-DBF/M);
    
    
    %create graph, use euclidean distance for edge weights
    %26-connectivity
    
    [nhood(:,1) nhood(:,2) nhood(:,3)] = ind2sub([3 3 3], 1:27);
    nhood = nhood - 2;
    nhood(all(nhood==0,2),:) = [];
    hood_weight = sqrt(sum(nhood.^2,2));
    hood_length = size(nhood,1);
    
    dest_node = zeros(N,hood_length);
    my_node = zeros(N,hood_length);
    edge_weight = zeros(N,hood_length);
    
    %substantially more efficient implementation of the below commented
    %area
    for n = 1:hood_length
        my_node(:,n) = 1:N;        
        
        dest_loc = node2loc + ones(N,1)*nhood(n,:);
        is_valid = all(dest_loc > 0, 2) & all(dest_loc <= ones(N,1)*im_max, 2);
        dest_ind = zeros(N,1);
        dest_ind(is_valid) = sub2ind(size(bin_im), dest_loc(is_valid,1), dest_loc(is_valid,2), dest_loc(is_valid,3));
        dest_node(is_valid, n) = ind2node(dest_ind(is_valid));
        is_valid = dest_node(:,n) ~= 0;
        edge_weight(is_valid, n) = hood_weight(n) * p_v(node2ind(dest_node(is_valid,n)));
    end
    %bottleneck block
%     for k = 1:N
%         my_node(k,:) = k;
% %         my_ind = node2ind(k);        
%         my_loc = node2loc(k,:);
% %         [my_loc(1) my_loc(2) my_loc(3)] = ind2sub(im_max, my_ind);            
%         
%         for n = 1:hood_length
%             dest_loc = my_loc + nhood(n,:);
%             if all(dest_loc > 0) && all(dest_loc <= im_max)
%                 
%                 dest_node(k,n) = loc2node(dest_loc(1), dest_loc(2), dest_loc(3));
%                     
%                 if dest_node(k,n) ~= 0
%                     edge_weight(k,n) = hood_weight(n) * p_v(node2ind(dest_node(k,n)));
%                     
%                     %edge_weight(k,n) = hood_weight(n) + p_v(node2ind(dest_node(k,n)));                    
%                     % ^ is as written in the paper. This makes little
%                     % sense, imo, since it just makes the actual length of
%                     % the path irrelevant
%                     
%                 end
%             end
%         end
%     end
%     toc
%     disp(sum(edge_weight(:)));


    has_edge = dest_node(:) ~= 0;
    G = sparse(my_node(has_edge), dest_node(has_edge), edge_weight(has_edge), N, N);
    
    
    %the commented code below is intended to find the soma and choose it as
    %the root node. For general purpose, we'll instead choose an arbitrary
    %point that's probably a leaf or very near one.
        
    
%     strat = zeros(im_max(1),1);
%     for c = 1:im_max(1)
%         strat(c) = sum(bin_im(c,:));
%     end
% 
%     cumstrat = cumsum(strat);
%     cumstrat = cumstrat/cumstrat(end);
%     
%     if ~exist('root_position','var') || isempty(root_position)    
%         %guess soma loc from stratification
%         
%         strat_med = find(cumstrat>.5,1,'first');
%         if strat_med > im_max(1)/2
%             root_position = 'GCL';
%         else
%             root_position = 'INL';
%         end
%         
%     end
%     
%     if strcmp(root_position, 'GCL');
%         slice_pos = find(cumstrat>=.05,1,'first');
%     else
%         slice_pos = find(cumstrat>=.95,1,'first');        
%     end
%         
%     im_slice = squeeze(bin_im(slice_pos,:,:));
%     [x_coord y_coord] = find(im_slice);
%     x_mean = mean(x_coord(:));
%     y_mean = mean(y_coord(:));
%     
%     dist_from_mean = (x_coord-x_mean).^2 + (y_coord-y_mean).^2;
%     [~, min_loc] = min(dist_from_mean);
%     
%     root_location = [slice_pos, x_coord(min_loc), y_coord(min_loc)];
%     root_node = ind2node(quick_ind(im_max, root_location));
    

    root_ind = find(bin_im(:),1,'first');
    root_node = ind2node(root_ind);
    
    root_nodes = root_node;
    
    while 1;        
        is_disconnected = true(1, size(G,1));
        for t = 1:length(root_nodes)
            is_disconnected = is_disconnected & isinf(graphshortestpath(G,root_nodes(t)));
        end
        
        if any(is_disconnected)
            root_nodes(end+1) = find(is_disconnected,1,'first');
        else
            break
        end
        
    end
        
    
    cube_out.nodes = [];
    cube_out.edges = [];
    cube_out.root_node = [];
    cube_out.node_radii = [];
    
    for t = 1:length(root_nodes)
        root_node = root_nodes(t);
    
        PDRF_G = graphshortestpath(G,root_node);
        PDRF = zeros(im_max);
        PDRF(node2ind) = PDRF_G;

        paths = [];
        k = 1;

        mask_im = bin_im & ~isinf(PDRF);

        while any(mask_im(:))

            mask_inds = find(mask_im(:));

            [~, farthest_in_ind] = max(PDRF(mask_inds));
            farthest_node = ind2node(mask_inds(farthest_in_ind));

            [~, paths{k}] = graphshortestpath(G,root_node,farthest_node);

            %prob a faster way to do this, but lazy and not a bottleneck
            for c = 1:length(paths{k})
                my_ind = node2ind(paths{k}(c));
                [my_loc(1) my_loc(2) my_loc(3)] = ind2sub(im_max, my_ind);

                r = DBF(my_ind)*scale + const;

                r_begin = ceil(my_loc - r);
                r_begin(r_begin<1) = 1;

                r_end = floor(my_loc + r);
                r_end = min([r_end; im_max]);

                mask_im(r_begin(1):r_end(1), r_begin(2):r_end(2), r_begin(3):r_end(3)) = false;
            end

            k = k + 1;

        end



        %consolidate paths

        tree_nodes = [];
        for k = 1:length(paths);
            tree_nodes = [tree_nodes; paths{k}'];
        end
        tree_nodes = unique(tree_nodes);
        tree_size = length(tree_nodes);

        inv_tree_nodes = sparse(tree_nodes, ones(tree_size,1), (1:tree_size)', N,1);
        root_node = inv_tree_nodes(paths{1}(1));

        edges = []; %probably losing a lot of speed here by not preallocating, but lazy
        for k = 1:length(paths);
            edges = [edges; full([inv_tree_nodes(paths{k}(1:end-1)') inv_tree_nodes(paths{k}(2:end)')])];
        end
        edges = unique(edges, 'rows');



        nodes = [];
        [nodes(:,1) nodes(:,2) nodes(:,3)] = ind2sub(im_max, node2ind(tree_nodes));
        node_radii = zeros(size(nodes,1),1, 'single');
        for n = 1:size(nodes,1)
            node_radii(n) = DBF(nodes(n,1), nodes(n,2), nodes(n,3));
        end

        cube_out.edges = [cube_out.edges; edges+size(cube_out.nodes,1)];
        cube_out.nodes = [cube_out.nodes; nodes];
        cube_out.root_node = [cube_out.root_node; root_node];
        cube_out.node_radii = [cube_out.node_radii; node_radii];
    
    end
end
        
function dup_nodes = find_duplicate_entries(nodes)    
    
    [sorted_nodes, sort_inds] = sortrows(nodes);
%     [dummy, unsort_inds] = sort(sort_inds);
    
    dup_nodes = zeros(size(nodes,1), 2);
    k = 1;
    for n = 1:size(sorted_nodes,1)-1
        if all(sorted_nodes(n,:) == sorted_nodes(n+1,:))
            dup_nodes(k,1:2) = [n, n+1];
            k = k+1;
        end
    end
    dup_nodes = dup_nodes(1:k-1,:);
    dup_nodes = sort_inds(dup_nodes);
    if size(dup_nodes,2) == 1
        dup_nodes = dup_nodes';
    end
    
    dup_nodes(dup_nodes(:,1) > dup_nodes(:,2),:) = dup_nodes(dup_nodes(:,1) > dup_nodes(:,2),[2 1]);
    
    dup_nodes = sortrows(dup_nodes);
end
    
    
    