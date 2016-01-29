function new_cns = subarbor_cell_creation(cell_num, guidebox, contact_threshold)

    C = get_constants;

    c_d = cell_data(cell_num);
    sp = c_d.get_surface;
    conns = double(c_d.contacts);
    soma_point = c_d.get_midpoint(true);
    
    new_cns = (1:9) * 10^floor(log(cell_num)/log(10)+1)+cell_num;
    for n = 1:9
        fh = figure; plot_cells(cell_num, 1, .01, [0 0 0]);
        plot(guidebox(:,1), guidebox(:,2), 'r');
        
        [x1, y1] = getpts(fh);
        
        if isempty(x1)
            break
        end
        
        x1 = [x1; x1(1)]; 
        y1 = [y1; y1(1)];
        
        in_p = inpolygon(sp(:,2), sp(:,3), x1, y1);
        
        surface_points = sp(in_p, :);
        p = zeros(100,3); 
        
        nearby_inds_cell = find_nearby_points(conns(3:5,:)', surface_points, contact_threshold);
        nearby_inds_isempty = false(size(nearby_inds_cell));
        for t = 1:length(nearby_inds_cell)
            nearby_inds_isempty(t) = isempty(nearby_inds_cell{t});
        end
        
        
        new_conns = conns(:,~nearby_inds_isempty);
        
        sfn = [C.point_dir, 'cell_' num2str(n), num2str(cell_num), '_surface.mat'];
        vfn = [C.raw_point_dir, 'points_', num2str(n), num2str(cell_num), '.mat'];
        cdfn = [C.cell_data_dir, 'cell_', num2str(n), num2str(cell_num), '_data.mat'];
        midfn = [C.soma_dir, 'cell_', num2str(n), num2str(cell_num), '_soma.mat'];
        
        save(sfn, 'surface_points')
        save(vfn, 'p');
        
        cd = cell_data(n*10^floor(log(cell_num)/log(10)+1)+cell_num);
        cd.contacts = new_conns;
        
        save(cdfn, 'cd');
        save(midfn, 'soma_point');
        close(fh)
    end
    new_cns = new_cns(1:n-1);
    
end
        