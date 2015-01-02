function write_points_binary(cell_num, fn)

    C = get_constants;

    f = fopen(fn,'w');
    
    point_files = get_files_with_names_including(C.raw_point_dir,num2str(cell_num));
            
    
    for n = 1:length(point_files)
        load([C.raw_point_dir point_files{n}]);
        fwrite(f, p', 'uint32');
    end
    fclose(f);        
end    
    