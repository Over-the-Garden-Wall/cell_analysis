function StratMap = make_stratification_map


    C = get_constants;
       
    
    point_dirs = dir(C.point_dir);    
    StratMap = containers.Map(1, zeros(100,1));
    
    for n = 1:length(point_dirs)
        if ~isempty(strfind(point_dirs(n).name,'.mat'))
            fn = [C.point_dir '/' point_dirs(n).name];
            
            uscore_locs = find(fn == '_');
            cell_num = str2double(fn(uscore_locs(end-1)+1 : uscore_locs(end) - 1));
            c_d = cell_data(cell_num);
            
            p = c_d.get_surface;
            d = C.f(p(:,1));
            strat = hist(d, -1:101);
            strat = strat(2:end-1);
            strat = strat/sum(strat);
            
            StratMap(cell_num) = strat;  
            disp(['Completed ' num2str(n) ' of ' num2str(length(point_dirs))]);
        end
       
    end
    
end