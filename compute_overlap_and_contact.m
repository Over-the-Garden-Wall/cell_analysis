function [overlap contact slope] = compute_overlap_and_contact(cell_group_A, cell_group_B)

    num_cells = [length(cell_group_A), length(cell_group_B)];
    
    overlap = zeros(num_cells);
    contact = zeros(num_cells);    
    
    for m = 1:num_cells(1)
        for n = 1:num_cells(2)
            
            overlap(m,n) = compute_cell_overlap([cell_group_A(m), cell_group_B(n)]);
            
            cell_dat = cell_data(cell_group_A(m));
            if cell_dat.contact_map.isKey(cell_group_B(n));
                contact(m,n) = cell_dat.contact_area(cell_dat.contact_map(cell_group_B(n)));
            else
                contact(m,n) = 0;
            end
            
            
        end
    end
    
    is_valid = overlap~=0;
    slope = mean(contact(is_valid)./overlap(is_valid));
    
end
    