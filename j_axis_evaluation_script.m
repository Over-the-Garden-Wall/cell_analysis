function j_axis_evaluation_script

    j_cell = 10010;
    
    C = get_constants;
    
    cell_types = {{'t1'}, {'t4'}, {'t3a', 't3b'}};
    
    
    
    axis_types = {'soma-distal', 'soma-pc', 'mean-distal', 'dist-distal', 'mean-pc'};
    
    
    j_dat = cell_data(j_cell);
    
    j_points = j_dat.get_surface;
    z_ax = C.f(j_points(:,1));
    j_points = j_points(z_ax>0 & z_ax<60, :);
    
    
    new_j_points = zeros(size(j_points));
    
    colors = colormap('Lines');
    
    for c = 1:length(cell_types)
        cell_nums = [];
        %get contact info
        for k = 1:length(cell_types{c})
            cell_nums = [cell_nums C.type.(cell_types{c}{k})];
        end
        
        
        num_cells = length(cell_nums);
        
        contacts = zeros(num_cells,1);
        locations = zeros(num_cells,3);
        
        line_h = zeros(length(axis_types),1);
        
        for n = 1:length(cell_nums)
            

            if j_dat.contact_map.isKey(cell_nums(n));
                cont_id = j_dat.contact_map(cell_nums(n));

                contacts(n) = j_dat.contact_count(cont_id);
            end

            cell_dat = cell_data(cell_nums(n));
            locations(n,:) = cell_dat.get_midpoint(false);

        end
        
        figure;
        gca; hold all;
        title(cell_types{c}{1});
        
        for n = 1:length(axis_types);
            
            switch axis_types{n}(1:4)
                case 'soma'
                    x0 = j_dat.get_midpoint(true);
                case 'mean'
                    x0 = j_dat.get_midpoint(false);
                case 'dist'
                    mp = j_dat.get_midpoint(true);  
                    x0 = mp;
                    old_dist = 0;
                    max_dist = -1;
                    while old_dist ~= max_dist
                        old_dist = max_dist;
                        x1 = x0;
                        dist = (x0(2)-j_points(:,2)).^2 + (x0(3)-j_points(:,3)).^2;
                        [max_dist max_ind] = max(dist);
                        x0 = j_points(max_ind,:);
                    end
                    
                    x0dist = (x0(2) - mp(2))^2 + (x0(3) - mp(3))^2;
                    x1dist = (x1(2) - mp(2))^2 + (x1(3) - mp(3))^2;
                    
                    if x1dist < x0dist
                        x0 = x1;
                    end
            end
            
        
            for d = 1:3
                new_j_points(:,d) = j_points(:,d)-x0(d);
            end        
            
            
            switch axis_types{n}(6:end)
                case 'distal'    
                    dist = new_j_points(:,2).^2 + new_j_points(:,3).^2;
                    [max_dist max_ind] = max(dist);
                    x1 = new_j_points(max_ind,:);                                        
                case 'pc'
                    
                    M = new_j_points'*new_j_points/size(j_points,1);
                    [eigvecs Ceigvals] = eigs(M);
                    x1 = eigvecs(:,1);
            end
            x1 = x1(2:3);
            x1 = x1/norm(x1);
            
            
            loc_on_axis = (locations(:,2)-x0(2)).*x1(1) + ...
                (locations(:,3)-x0(3)).*x1(2);                        
                        
            loc_on_axis = loc_on_axis/1000;
            
            [loc_on_axis sort_ind] = sort(loc_on_axis);
            contacts = contacts(sort_ind);
            
            line_h(n) = scatter(loc_on_axis, contacts, '*', 'MarkerEdgeColor',colors(n,:),'MarkerFaceColor',colors(n,:));
            
            
            linearCoef = polyfit(loc_on_axis,contacts,4);
            linearFit = polyval(linearCoef,loc_on_axis);
            plot(loc_on_axis,linearFit,'-', 'Color', colors(n,:));
                        
        end
        legend(line_h, axis_types);
    end
end
                    