function [mean_contact, ste_contact] = quick_contact_survey(ref_cells, cont_cells, ref_range)

    if ~exist('ref_range', 'var') || isempty(ref_range)
        ref_range = [-Inf Inf];
    end
    
    num_ref = length(ref_cells);
    num_cont = length(cont_cells);
    mean_contact = zeros(num_cont, 1);
    ste_contact = zeros(num_cont ,1);
    
    ref_mids = zeros(num_ref,2);
    for r = 1:num_ref
        c_d = cell_data(ref_cells(r));
        m = c_d.get_midpoint;
        ref_mids(r,:) = m(2:3);
    end
    
    for c = 1:num_cont
        c_d = cell_data(cont_cells(c));
        
        my_contacts = double(c_d.contacts);
        
        is_valid = false(1, size(my_contacts,2));
        for r = 1:num_ref
            is_me = find(my_contacts(1,:)==ref_cells(r));
            sub_conts = my_contacts(:,is_me);
            d = sqrt((sub_conts(4,:) - ref_mids(r,1)).^2 + (sub_conts(5,:) - ref_mids(r,2)).^2);
            is_valid(is_me(d>=ref_range(1) & d<=ref_range(2))) = true;
        end
        
        mean_contact(c) = mean(my_contacts(2,is_valid));
        ste_contact(c) = std(my_contacts(2,is_valid)) / sqrt(sum(is_valid)-1);
        
    end
    
    
end
        