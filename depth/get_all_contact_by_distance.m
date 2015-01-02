function [dist, cont_area, soma_point] = get_all_contact_by_distance(target_cell, contacting_cells, conns, dist_mode)



    C = get_constants;
    num_cells = length(contacting_cells);
    
    cont_area = cell(num_cells,1);    
    dist = cell(num_cells,1);
    
    fn = [C.point_dir '/cell_' num2str(target_cell) '_surface.mat'];
    load(fn);
    
    soma_point = get_soma_loc(target_cell); 
    
    if ~exist('dist_mode','var') || ~strcmp(dist_mode, 'axis')
        warning('unknown dist_mode, defaulting to ''point''');
    elseif strcmp(dist_mode, 'axis')        
        conns(4,:) = 0;
        soma_point(1) = 0;    
    end
    
    sub_conns = pick_conns(conns, target_cell, contacting_cells);
    
        
    for k = 1:num_cells
        num_contacts = size(sub_conns{k},2);
        dist{k} = sqrt(sum((double(sub_conns{k}(4:6,:)) - soma_point'*ones(1,num_contacts)).^2,1));            
    end                

    
    for k = 1:num_cells; 
        cont_area{k} = double(sub_conns{k}(3,:));
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