function section_lengths = get_lengths_for_cell_blocks(cell_no, cell_name)
    dsmp_fact = [3 3 3];
    im = points2full(cell_no, dsmp_fact);
    
    [nodes edges] = skeletonize_by_grid(im,[32 32 32]);
    [nodes edges] = condense_skeleton(nodes, edges, 5);
    
    
    
    load(['./forSrini/' cell_name '_info.mat']);
    
    section_lengths = zeros(length(omni_file_fns),1);
    nodes = nodes*3;
    
    for n = 1:length(omni_file_fns)
        [start_coords, end_coords] = omni_fn2coords(omni_file_fns{n});
        [subnodes, subedges] = prune_skeleton_to_block(nodes,edges,start_coords, end_coords);
        
        section_lengths(n) = calc_skele_length(subnodes, subedges);
    end
end