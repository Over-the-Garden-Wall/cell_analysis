function [total_contact, total_vox_in_hull] = get_contact_density_whulls(contact_cells, hull_cells)


    C = get_constants;
    
    num_hc = length(hull_cells);
    num_cc = length(contact_cells);
    
    hulls = cell(num_hc,1);
    for n = 1:num_hc
        cell_dat = cell_data(hull_cells(n));        
        hulls{n} = cell_dat.hull_2d;
    end
    
    
    total_contact = zeros(num_cc,num_hc);
    total_vox_in_hull = zeros(num_cc,num_hc);
    
    for c = 1:num_cc
        cell_dat = cell_data(contact_cells(c));
        p = cell_dat.get_surface;
        p = p(:,2:3);
        
        for h = 1:num_hc
                if cell_dat.contact_map.isKey(hull_cells(h))
                    map_id = cell_dat.contact_map(hull_cells(h));
                    total_contact(c,h) = cell_dat.contact_area(map_id);
                end
                
                if ~exist([C.hull_dir 'c' num2str(contact_cells(c)) '/'], 'dir')
                    mkdir([C.hull_dir 'c' num2str(contact_cells(c)) '/']);
                end
                
                hull_fn = [C.hull_dir 'c' num2str(contact_cells(c)) '/c' num2str(hull_cells(h)) '.mat'];
                
                
            if exist(hull_fn, 'file')
                
                load(hull_fn);                
                
                
            else
            

                in_p = inpolygon(p(:,1),p(:,2),hulls{h}(:,1),hulls{h}(:,2));

                points_in_hull = sum(in_p);
                save(hull_fn, 'points_in_hull');
            end
            total_vox_in_hull(c,h) = points_in_hull;
                
                
            
        end
    end
end
                
            
            
            