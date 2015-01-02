function [mean_loc total_contact cont_nums] = get_axial_contacts(cell_num, contact_cells, varargin)
    
    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addRequired('contact_cells', @isnumeric);
    p.addOptional('trunk_cutoff', [-Inf Inf], @isnumeric);
    p.addOptional('use_cont_midpoint', false, @islogical);
    
    
    p.parse(cell_num, contact_cells, varargin{:});    
    s = p.Results;
    
    C = get_constants;
    
    
    ref_cell = cell_data(cell_num);
    
    ref_mean = ref_cell.get_midpoint(true);
    
    num_conts = length(contact_cells);
    
    
    contacts = double(ref_cell.contacts);
    
    is_valid = false(1,size(contacts,2));
    for k = 1:num_conts
        is_valid = is_valid | contacts(1,:) == contact_cells(k);
    end
    
    total_contact = contacts(2, is_valid)';
    cont_nums = contacts(1,is_valid);        
    
    cont_loc = contacts(3:5,is_valid)';

    if ~isempty(cont_loc)
        is_valid = C.f(cont_loc(:,1)) >= s.trunk_cutoff(1) & C.f(cont_loc(:,1)) <= s.trunk_cutoff(2);
        total_contact = total_contact(is_valid);
        cont_loc = cont_loc(is_valid,:);
        cont_nums = cont_nums(is_valid);
    end
    if ~s.use_cont_midpoint
        cont_loc = zeros(length(cont_nums),3);
        for n = 1:length(cont_nums)
            cell_dat = cell_data(cont_nums(n));
            cont_loc(n,:) = cell_dat.get_midpoint(false);            
        end
    end
    
    
    for d = 1:3
        cont_loc(:,d) = cont_loc(:,d) - ref_mean(d);
    end
    
    
    
            if ~ref_cell.is_symmetric
                mean_loc = cont_loc(:,2)*ref_cell.dist_axis(1) + cont_loc(:,3)*ref_cell.dist_axis(2);
            else
                mean_loc = sqrt(cont_loc(:,2).^2+cont_loc(:,3).^2);
            end
    
end