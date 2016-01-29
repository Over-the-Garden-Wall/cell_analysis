function h = visualize_contacts(ref_cell, sac_range, contact_minimum)

    %code incomplete
    
    
    if ~exist('sac_range', 'var') || isempty(sac_range)
        sac_range = [0 Inf];
    end
    if ~exist('contact_minimum', 'var') || isempty(contact_minimum)
        contact_minimum = 100;
    end

    C = get_constants;
    
    figure; hold all
    plot_cells(ref_cell, 1, .01, [0 0 0]);
    
    
    c_d = cell_data(ref_cell);
    conts = double(c_d.contacts);
    for n = 1:length(cell_nums)
        
        is_me = false(1,size(conts,2));
        for k = 1:size(conts,2)
            if any(conts(1,k)==cell_nums{n})
                is_me(k) = true;
            end
        end
        h(n) = scatter(conts(4,is_me), conts(5,is_me), '*');
    end
    