function visualize_conn_locs(cell_num, cont_nums, cont_threshold)

    if ~iscell(cont_nums)
        cont_nums = {cont_nums};
    end
    
    cell_dat = cell_data(cell_num);
    
    
    figure; hold all;
    plot_cells(cell_num,1,.01,[0 0 0]);
    
    cont_locs = cell(length(cont_nums),1);    
    for t = 1:length(cont_nums)
        for c = 1:length(cont_nums{t})
            cn = cont_nums{t}(c);
            sub_conn = double(cell_dat.contacts(:,cell_dat.contacts(1,:)==cn & cell_dat.contacts(2,:)>=cont_threshold));
            cont_locs{t} = [cont_locs{t}; sub_conn(4:5,:)'];
        end
        
        scatter(cont_locs{t}(:,1), cont_locs{t}(:,2), '*');
    end
    
    
end