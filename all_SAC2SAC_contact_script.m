function density = all_SAC2SAC_contact_script

    C = get_constants;
    
    num_sacs = length(C.type.off_sac);
    
    density = cell(num_sacs,1);
    
    size_max = [0 0 0];
    for n = 1:num_sacs
        d{n} = get_contacts_all(C.type.off_sac(n), C.type.off_sac, true);
        size_max = max([size_max; size(d{n})]);
    end
    
    density = zeros(size_max);
    for n = 1:num_sacs
        dSz = size(d{n});
        density(1:dSz(1), 1:dSz(2), 1:dSz(3)) = density(1:dSz(1), 1:dSz(2), 1:dSz(3)) + ...
            d{n};
    end
    
end
    
    