function visualize_sac_connectivity(sac_num)

    C = get_constants;

    sacs = [C.type.on_sac];%, C.type.off_sac];

    
    c_d = cell_data(sac_num);
    contacts = double(c_d.contacts);
    
    is_sac = false(max(contacts(1,:)),1);
    is_sac(sacs) = true;
    
    contacts = contacts(:,is_sac(contacts(1,:)));
    
    contact_sacs = unique(contacts(1,:));
    contact_total = zeros(size(contact_sacs));
    
    for s = 1:length(contact_sacs);
        contact_total(s) = sum(contacts(2,contacts(1,:)==contact_sacs(s)));
    end
%     contact_total = contact_total/max(contact_total);
    contact_total = contact_total > median(contact_total);
    
    figure; plot_cells(sac_num, 1, .01, ones(1,3)*.5);
    
    for s = 1:length(contact_sacs)
        c_d = cell_data(contact_sacs(s));
        m = c_d.get_midpoint(true);
        clr = [1 0 0] * contact_total(s) + (1-contact_total(s)) * [0 0 1];
        
        scatter(m(2), m(3), '*', 'markerEdgeColor', clr);
    end
    
    
end