function id_table = transform_and_plot(cell_nums, Q, color_pick)
    res = [16.5 16.5 25];
    num_cells = length(cell_nums);
    
    y_data = cell(num_cells,1);
    x_data = cell(num_cells,1);
    id_table = zeros(num_cells, 6);
    id_table(:,1) = cell_nums;
    
%     f = inline(' (30-62)/(2.7-1.65)/10^4.*(x-2.7*10^4)+30');
    f = inline(' (30-108)/(2.7-1.5)/10^4.*(x-2.7*10^4)+30');
    
    
    for n = 1:num_cells
                
        
        fn = ['./surface_points/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            for k = 1:3
                surface_points(:,k) = surface_points(:,k)*res(k);
            end                
            
            surface_points = surface_points*Q';
            
            [y_data{n} x_data{n}] = hist(surface_points(:,1),20);
            y_data{n} = y_data{n}/sum(y_data{n});
            
            [id_table(n,3), max_ind] = max(y_data{n});
            color_val = x_data{n}(max_ind);
            
            id_table(n,4:6) = [1 0 0]*(1-abs(20000-color_val)/10000) + ...
                [0 1 0]*(1-abs(25000-color_val)/5000) + ...
                [0 0 1]*(1-abs(30000-color_val)/10000);
            id_table(n,4:6) = max([id_table(n,4:6); 0 0 0]);
            id_table(n,4:6) = min([id_table(n,4:6); 1 1 1]);
            
            
            %change x axis!
            x_data{n} = f(x_data{n});
            id_table(n,2) = x_data{n}(max_ind);            
            
                                    
        end
    end
    
    figure;        
    is_valid = true(num_cells,1);
    for n = 1:num_cells

        [dummy b] = max(id_table(n,4:6));
        if (color_pick == 0 || color_pick == b) && ~isempty(x_data{n})


            plot(x_data{n}, y_data{n}, 'Color', id_table(n,4:6));
            hold on;
        else
            is_valid(n) = false;
        end

    end
    
    figure;     
    valid_list = find(is_valid');
    p = ceil(sqrt(length(valid_list)));
    q = ceil(length(valid_list)/p);
    
    for k = 1:length(valid_list)
        subplot(p,q,k); 
        n = valid_list(k);
        
        
        plot(x_data{n}, y_data{n}, 'Color', id_table(n,4:6));
        
        title(num2str(cell_nums(n)))
%             hold on;

    end
    
    
    id_table(~is_valid,:) = [];
    [dummy sort_inds] = sort(id_table(:,2));
    id_table = id_table(sort_inds,:);
    
    disp({'cell_id', 'ipl_depth', 'max_height', 'color'}) 
    for n = 1:size(id_table)
        disp([num2str(id_table(n,1)) '   ' num2str(id_table(n,2)) '   ' num2str(id_table(n,3)), '   ' ...
            '[' num2str(id_table(n,4)) ', ' num2str(id_table(n,5)) ', ' num2str(id_table(n,6)) ']']);
    end
     
    
end