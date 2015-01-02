function contact_data = summarize_dendrites(sac_nums, path_nums, types)

    C = get_constants;

    allowed_hull = find_covered_area(types);

    contact_data = cell(length(sac_nums),1);

    for n = 1:length(sac_nums)
        
        contact_list = analyze_dendrite(sac_nums(n), path_nums(n), allowed_hull, 1000);
        
        cont_types = zeros(size(contact_list,1),1);
        for k = 1:size(contact_list,1);
            for t = 1:length(types)
                if any(C.type.(types{t})==contact_list(k,1))
                    cont_types(k) = t;
                end
            end
        end
        
        contact_data{n} = [contact_list(cont_types>0,:) cont_types(cont_types>0)];
        
    end
end

    