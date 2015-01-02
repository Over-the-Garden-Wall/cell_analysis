types_to_write = {'t1','t2','t3a','t3b','t4', 'sure_off_sac'};

C = get_constants;

for t = 1:length(types_to_write)
    for n = C.type.(types_to_write{t});
        write_points_binary(n,['~/stratification/point_data_raws/c' num2str(n) '.raw']);
    end
end