function [total_contact, total_vox_in_hull] = get_contact_density_whulls(contact_cells, hull_cells, depth_range, contact_range)

    if ~exist('depth_range','var') || isempty(depth_range)
        depth_range = [0 100];
    end
    if ~exist('contact_range','var') || isempty(contact_range)
        contact_range = [0 Inf];
    end


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
        d = C.f(p(:,1));
        p = p(d >= depth_range(1) & d <= depth_range(2),2:3);
        
        conts = double(cell_dat.contacts);
        conts = conts(:,conts(2,:) >= contact_range(1) & conts(2,:) <= contact_range(2));
        
        cont_depth = C.f(conts(3,:));
        
        for h = 1:num_hc
                
                total_contact(c,h) = sum(conts(2, conts(1,:) == hull_cells(h) & ...
                    cont_depth >= depth_range(1) & cont_depth <= depth_range(2)));
                
                
                cell_dir = [C.hull_dir 'c' num2str(contact_cells(c)) '/'];
                if ~exist(cell_dir, 'dir')
                    mkdir(cell_dir);
                end
                cell_hull_dir = [cell_dir 'c' num2str(hull_cells(h)) '/'];
                if ~exist(cell_hull_dir, 'dir')
                    mkdir(cell_hull_dir);
                end
                
                hull_fn = [cell_hull_dir 'd_' num2str(depth_range(1)), '_', num2str(depth_range(2)), '.mat'];
                
                
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
                
            
            
            