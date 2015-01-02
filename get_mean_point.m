function mean_point = get_mean_point(cell_num, use_soma)

    if use_soma
        mean_point = get_soma_loc(cell_num); 
    else
        mean_point = get_mean_loc(cell_num);
    end
    
end