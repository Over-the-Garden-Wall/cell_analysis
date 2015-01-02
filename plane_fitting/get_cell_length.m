function full_len = get_cell_length(cell_name)

    [nodes edges] = get_skeleton(cell_name);
    full_len = get_skele_length(nodes,edges);
    


    
    
    
end