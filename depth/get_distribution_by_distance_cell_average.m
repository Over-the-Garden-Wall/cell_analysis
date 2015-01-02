function [mean_dist, total_area, soma_point] = get_distribution_by_distance_cell_average(target_cell, contacting_cells, conns, dist_mode, loc_mode)



    C = get_constants;
    num_cells = length(contacting_cells);

    fn = [C.point_dir '/cell_' num2str(target_cell) '_surface.mat'];
    load(fn);
    
    soma_point = get_soma_loc(surface_points); %#ok<NODEF>
    
    if ~exist('dist_mode','var') || ~strcmp(dist_mode, 'axis')
        warning('unknown dist_mode, defaulting to ''point''');
    elseif strcmp(dist_mode, 'axis')        
        conns(4,:) = 0;
        soma_point(1) = 0;    
    end
    
    sub_conns = pick_conns(conns, target_cell, contacting_cells);
    
    mean_dist = zeros(num_cells,1);
    if ~exist('loc_mode','var') || ~strcmp(loc_mode, 'arbor')
        if ~strcmp(loc_mode, 'contact')
            warning('unknown dist_mode, defaulting to ''contact'''); %#ok<*WNTAG>
        end
        for k = 1:num_cells
            weights = double(sub_conns{k}(3,:));
            weights = weights/sum(weights);
            num_contacts = size(sub_conns{k},2);
            dist = sqrt(sum((double(sub_conns{k}(4:6,:)) - soma_point'*ones(1,num_contacts)).^2,2));
            mean_dist(k) = sum(dist.*weights);
        end
        
    elseif strcmp(dist_mode, 'arbor')
        for k = 1:num_cells

            fn = [C.point_dir '/cell_' num2str(target_cell) '_surface.mat'];
            load(fn);
        
            if strcmp(dist_mode, 'axis')        
                surface_points(:,1) = 0;
            else
                surface_points(:,1) = f(surface_points(:,1));            
            end
            
            num_points = size(surface_points,1);
            dist = sqrt(sum((surface_points - ones(num_points,1)*soma_point).^2));
            mean_dist(k) = mean(dist);
        end
        
    end
    
    

    total_area = zeros(num_cells,1);
    for k = 1:num_cells; 
        total_area(k) = double(sum(sub_conns{k}(3,:)));
    end
        
% 
% cols = [ones(0,1)*[1 0 0]; ones(100,1)*[0 0 1]];
% figure; hold on;
% for n = 1:num_cells;
%     scatter(mean_dist(n), total_area(n), '*', 'MarkerEdgeColor', cols(n,:), 'MarkerFaceColor', cols(n,:));    
% end
% cell_names = cell(num_cells,1);
% for n = 1:num_cells
%     cell_names{n} = num2str(contacting_cells(n));
% end
% 
% 
% 
% % figure;
% % scatter(mean_dist, total_area,'*');    
% text(mean_dist,total_area + max(total_area)/30, cell_names);
% % scatter(mean_dist, total_area);    
% 
% % x = x -;
% 
% % f = inline('-.006*(x-2.7*10^4)+30');
% % x = f(x);
% c = corr(mean_dist, total_area);
% 
% %     
% figure; hold on;
% max_val = 0;
% for n = 1:num_cells;
%     fn =['./surface_points/cell_' num2str(contacting_cells(n)) '_surface.mat'];
%     if exist(fn,'file')
%         point_info = whos('-file', fn);
% 
%         scatter(mean_dist(n), total_area(n)/point_info.size(1), 'MarkerEdgeColor', cols(n,:), 'MarkerFaceColor', cols(n,:));
%         max_val = max(max_val, total_area(n)/point_info.size(1));        
%     end
% end
% 
% for n = 1:num_cells;
%     fn =['./surface_points/cell_' num2str(contacting_cells(n)) '_surface.mat'];
%     if exist(fn,'file')
%         point_info = whos('-file', fn);
% 
%         text(mean_dist(n), total_area(n)/point_info.size(1) + max_val/30, cell_names{n}) %, 'MarkerEdgeColor', cols(n,:), 'MarkerFaceColor', cols(n,:));
%     end
% end


                
%     total_points = total_points + ;


end