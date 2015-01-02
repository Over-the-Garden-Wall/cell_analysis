function [mean_loc total_contact is_valid] = get_mean_axial_contacts(cell_num, contact_cells, varargin)
    
    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addRequired('contact_cells', @isnumeric);
    p.addOptional('use_count', false, @islogical);
    p.addOptional('use_soma', false, @islogical);
    p.addOptional('use_contacts', false, @islogical);
    p.addOptional('include_zeros', false, @islogical);
    
    
    p.parse(cell_num, contact_cells, varargin{:});    
    s = p.Results;
    
    C = get_constants;
    
    
    ref_cell = cell_data(cell_num);
    
    ref_mean = ref_cell.get_midpoint(true);
    
    num_conts = length(contact_cells);
    
    
    mean_loc = zeros(num_conts,1);
    total_contact = zeros(num_conts,1);
    is_valid = true(num_conts,1);    
    
    
    for k = 1:num_conts
        if ref_cell.contact_map.isKey(contact_cells(k)) || s.include_zeros;
            if ref_cell.contact_map.isKey(contact_cells(k))
                cont_num = ref_cell.contact_map(contact_cells(k));
                if s.use_count
                    total_contact(k) = ref_cell.contact_count(cont_num);
                else
                    total_contact(k) = ref_cell.contact_area(cont_num);
                end
            end
            
            off_cell = cell_data(contact_cells(k));
                       
            if s.use_contacts
                c = double(ref_cell.get_contacts);
                c = c(:,c(1,:)==contact_cells(k));
                off_loc = [sum(c(3,:).*c(2,:)), sum(c(4,:).*c(2,:)), sum(c(5,:).*c(2,:))];
                off_loc = off_loc/sum(c(2,:));
            else
                off_loc = off_cell.get_midpoint(s.use_soma);
            end
            
            cell_diff = off_loc-ref_mean;
            
            if ~ref_cell.is_symmetric
                mean_loc(k) = cell_diff(2)*ref_cell.dist_axis(1) + cell_diff(3)*ref_cell.dist_axis(2);
            else
                mean_loc(k) = sqrt(cell_diff(2)^2+cell_diff(3)^2);
            end
            
            
        else
            is_valid(k) = false;
        end
    end
    
end