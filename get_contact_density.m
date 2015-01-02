function [] = get_contact_density(contact_cells, hull_cells)


    C = get_constants;
    
    num_hc = length(hull_cells);
    num_cc = length(contact_cells);
    
    hulls = cell(num_hc,1);
    for n = 1:num_hc
        cell_dat = cell_data(hull_cells(n));        
        hulls{n} = cell_dat.get_2d_hull;
    end
    
    
    total_contact = zeros(num_cc,num_hc
    
    for c = 1:num_cc
        for h = 1:num_hc