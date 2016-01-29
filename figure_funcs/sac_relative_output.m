function [observed_contacts, expected_contacts] = sac_relative_output(cell_nums)
    
    C = get_constants;

    lat_res = 100;
    
    sac_cont_lateral = zeros(lat_res);
    sac_cont_vertical = zeros(100,1);
    
    
    sacs = [C.type.off_sac];% C.type.on_sac];
    sac_conts = cell(length(sacs),1);
    lat_max = [0 0];
    
    for s = 1:length(sacs);
        sac_conts{s} = detect_vericose_contacts(sacs(s), 500, 200, 50000);
        lat_max = max([sac_conts{s}(4:5,:)'; lat_max]);
    end
    
    for s = 1:length(sacs);
        d = round(C.f(sac_conts{s}(3,:)));
        d = d(d >= 1 & d <= 100);
        for n = d
            sac_cont_vertical(n) = sac_cont_vertical(n) + 1;
        end
        
        p = [ceil(sac_conts{s}(4,:) / lat_max(1) * lat_res); ...
            ceil(sac_conts{s}(5,:) / lat_max(2) * lat_res)]';
        for n = 1:size(p,1)
            sac_cont_lateral(p(n,1), p(n,2)) = sac_cont_lateral(p(n,1), p(n,2)) + 1;
        end
    end
    sac_cont_vertical = sac_cont_vertical / sum(sac_cont_vertical);
    sac_cont_vertical(sac_cont_vertical < .0001) = 0;
    
    num_cells = length(cell_nums);
    expected_contacts = zeros(num_cells,1);
    observed_contacts = zeros(num_cells,1);
    
    for cn = 1:num_cells
        
        try
            
        c = cell_nums(cn);
        
        c_d = cell_data(c);
        p = c_d.get_surface;
        
        d = round(C.f(p(:,1)));
        
        p(d<1 | d>100,:) = [];
        d(d<1 | d>100,:) = [];
        
        is_valid = sac_cont_vertical(d) > 0;
        p = p(is_valid,:);
        d = d(is_valid);
        
        p = ceil([p(:,2) / lat_max(1) * lat_res, p(:,3) / lat_max(2) * lat_res]);
        p(p>lat_res) = lat_res;
        
        for n = 1:length(d)
            expected_contacts(cn) = expected_contacts(cn) + sac_cont_vertical(d(n)) * sac_cont_lateral(p(n,1), p(n,2)) / 5000;
            %5000 is the oodsgc normalization factor, found empirically
        end
        
        for s = 1:length(sacs);
            observed_contacts(cn) = observed_contacts(cn) + sum(sac_conts{s}(1,:) == c);
        end
        
        
        catch ME
            disp(ME.message)        
        end
        
    end
    
    
    
end