function type_name = cell_type(cell_num)
    C = get_constants;

    fns = fieldnames(C.type);
    type_name = 'unknown';
    for n = 1:length(fns)
        if any(C.type.(fns{n}) == cell_num);
            type_name = fns{n};
            break
        end
    end
end