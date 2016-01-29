function make_stuff_for_jinseop


    C = get_constants;
    
    point_dirs = dir(C.point_dir);        
    cell_nums = zeros(length(point_dirs),1);
    for n = 1:length(point_dirs)
        if ~isempty(strfind(point_dirs(n).name,'.mat'))
            us = find(point_dirs(n).name == '_');
            cell_nums(n) = str2double(point_dirs(n).name(us(1)+1:us(2)-1));
        end
    end
    
    cell_nums(cell_nums == 0) = [];
    
    num_cells = length(cell_nums);
    
    cell_map = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for k = 1:length(cell_nums)
        cell_map(cell_nums(k)) = k;
    end
    
    save('/data/home/greenem/data/stratification/jinseop_helper/cell_map.mat','cell_map');
    
    load(C.conn_loc);
    conns = double(conns);
    
    save('/data/home/greenem/data/stratification/jinseop_helper/cell_contacts.mat','conns');
    
    all_strats = zeros(100, num_cells);
    for k = 1:length(cell_nums)
        c_d = cell_data(cell_nums(k));
        s = c_d.stratification;
        s(1:20) = [];
        if length(s) > 100
            s = s(1:100);
        end
        
        all_strats(1:length(s),k) = s/sum(s);
    end
                
    save('/data/home/greenem/data/stratification/jinseop_helper/cell_strats.mat','all_strats');
    
    cell_hulls = cell(num_cells,1);
    for k = 1:num_cells
        c_d = cell_data(cell_nums(k));
        [cell_hulls{k}(:,1), cell_hulls{k}(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
    end
        
    save('/data/home/greenem/data/stratification/jinseop_helper/cell_hulls.mat','cell_hulls');
end
    
    