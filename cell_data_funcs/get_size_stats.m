function [SA, vol] = get_size_stats(cell_nums)

    C = get_constants;

    vol_dir = C.raw_point_dir;
    
    sa_dir = C.point_dir;

    num_cells = length(cell_nums);
    SA = zeros(num_cells,1);
    vol = zeros(num_cells,1);
    
    for n = 1:num_cells

        fn =[sa_dir '/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            point_info = whos('-file', fn);
            SA(n) = point_info.size(1);

            fns = get_files_with_names_including(vol_dir, num2str(cell_nums(n)));

            for k = 1:length(fns)
                point_info = whos('-file', [vol_dir fns{k}]);
                vol(n) = vol(n)+point_info.size(1);
            end



        end
    
    end
end