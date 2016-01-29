function sac_to_sac_central_contacts

    C = get_constants;
    size_threshold = 50;
    min_distance_from_boundary = 150000;
    num_angle_bins = 3;
    distal_threshold = 66000;
    
%     silly_dist_threshold = 25000;
    
    angle_bins = (.5:num_angle_bins-.5) / num_angle_bins * pi;
    
    bounds = [C.min_xy; C.max_xy];
    valid_region = [bounds(1,:) + min_distance_from_boundary; ...
        bounds(2,:) - min_distance_from_boundary];
    
    l_names = {'sure_off_sac', 'on_sac'};
    
    list_names = {'cell_ids', 'contact_loc', 'contact_size'};
    
%     cn_set = {[70014, 70158, 70016, 70154], [70164, 70176, 70230, 70205], ...
%         [70130, 70080, 70147, 70116], [70197, 70169, 70244, 70191], ...
%         [70102 70134 70027 70089], [70211 70215 70216 70188]};
    cn_set = {C.type.sure_off_sac, C.type.on_sac};
%     cn_set = {70014, 70158, 70016, 70154}; 
%     , [70164, 70176, 70230, 70205], ...
%     
% [70130, 70080, 70147, 70116],
% [70197, 70169, 70244, 70191], ...
%             [70102 70134 70027 70089], 
%         [70211 70215 70216 70188]};

%     for l = 1:2
    for l = 1:length(cn_set)
        
%     cns =  [70014, 70158, 70016, 70154]; %off center
%     cns =  [70164, 70176, 70230, 70205]; %on center
% cns = [70130, 70080, 70147, 70116]; %off l-edge
% cns = [70197, 70169, 70244, 70191]; %on l-edge
% cns = [70102 70134 70027 70089]; %off bl-corner
% cns = [70211 70215 70216 70188]; %on bl-corner

        
%         l_info{l} = contact_info(C.type.(l_names{l}), C.type.(l_names{l}), true);
%         l_info{l} = contact_info(70014, C.type.(l_names{l}));
%         l_info{l} = contact_info(cns, C.type.(l_names{l}));
        l_info{l} = contact_info(cn_set{l}, [C.type.(l_names{1}), C.type.(l_names{2})], true);

        valid_contacts = l_info{l}.cell_ids(:,1) ~= l_info{l}.cell_ids(:,2);
        valid_contacts = valid_contacts & l_info{l}.contact_size > size_threshold;
        valid_contacts = valid_contacts & ...
            (l_info{l}.contact_loc(:,2) > valid_region(1,1) & ...
            l_info{l}.contact_loc(:,3) > valid_region(1,2) & ...
            l_info{l}.contact_loc(:,2) < valid_region(2,1) & ...
            l_info{l}.contact_loc(:,3) < valid_region(2,2));
%         valid_contacts = valid_contacts & valid_somata(l_info{l}.cell_ids(:,1)) & ...
%             valid_somata(l_info{l}.cell_ids(:,2));
                                                               
        for k = 1:length(list_names);
            l_info{l}.(list_names{k}) = l_info{l}.(list_names{k})(valid_contacts,:);
        end
        
        soma_locs = [l_info{l}.soma_loc(l_info{l}.cell_ids(:,1),2:3), ...
            l_info{l}.soma_loc(l_info{l}.cell_ids(:,2),2:3)];
        
        rel_locs = [l_info{l}.contact_loc(:,2:3), l_info{l}.contact_loc(:,2:3)] - soma_locs;
        dists = [sqrt(sum(rel_locs(:,1:2).^2,2)), sqrt(sum(rel_locs(:,3:4).^2,2))];
        
        rel_locs(:,1) = rel_locs(:,1)./dists(:,1);
        rel_locs(:,2) = rel_locs(:,2)./dists(:,1);
        rel_locs(:,3) = rel_locs(:,3)./dists(:,2);
        rel_locs(:,4) = rel_locs(:,4)./dists(:,2);
        angles = acos(rel_locs(:,1).*rel_locs(:,3) + rel_locs(:,2).*rel_locs(:,4));   
        
%         angles(abs(dists(:,1)-dists(:,2))<silly_dist_threshold) = [];
        
        layer_plot_data{l} = hist(angles, angle_bins);
%         layer_plot_data{l} = weighted_hist(angles, l_info{l}.contact_size, angle_bins);
        
        
        
%         dists = [sqrt(sum(rel_locs(:,1:2).^2,2)), sqrt(sum(rel_locs(:,3:4).^2,2))];
%         is_proximal = dists < distal_threshold;
%         is_pd = sum(is_proximal,2) == 1;
%         is_dd = sum(is_proximal,2) == 0;
%         
%         layer_dist_plot{l}(:,1) = hist(angles(is_pd), angle_bins);
%         layer_dist_plot{l}(:,1) = layer_dist_plot{l}(:,1)/sum(layer_dist_plot{l}(:,1));
%         layer_dist_plot{l}(:,2) = hist(angles(is_dd), angle_bins);
%         layer_dist_plot{l}(:,2) = layer_dist_plot{l}(:,2)/sum(layer_dist_plot{l}(:,2));
%         layer_dist_plot{l}(:,3) = hist(angles(~is_dd & ~is_pd), angle_bins);
%         layer_dist_plot{l}(:,3) = layer_dist_plot{l}(:,3)/sum(layer_dist_plot{l}(:,3));
        
    end
    
%     linespecs = {'-o', '--o', '-s', '--s', '-^', '--^'};
    linespecs = {'-', '-'};
%     linespecs = {'-', '-', '-', '-'};

    figure; hold all
    for l = 1:length(layer_plot_data)
%         plot(angle_bins, layer_plot_data{l}/sum(layer_plot_data{l}), linespecs{l}, 'lineWidth', 2);
        plot(angle_bins, layer_plot_data{l}, linespecs{l}, 'lineWidth', 2);
    end
%     plot(angle_bins, layer_plot_data{2}/sum(layer_plot_data{2}), 'lineWidth', 2);
    
%     figure; hold all
%     plot(angle_bins, layer_plot_data{1}, 'lineWidth', 2);
%     
%     
%     figure; hold all
%     plot(angle_bins, layer_dist_plot{1}, 'lineWidth', 2);
%     plot(angle_bins, layer_dist_plot{2}, 'lineWidth', 2);
% %    
    

    dist_thresh = 80000;
    
    
    
%     valid_region = valid_region - 50000;
    %random points control
%     for l = 1:2
% %         num_somata = ceil(sqrt(size(l_info{l}.soma_loc,1)))^2;
% %         [x, y] = meshgrid((.5:sqrt(num_somata))/sqrt(num_somata)*(bounds(2,1)-bounds(1,1)) + bounds(1,1), ...
% %             (.5:sqrt(num_somata))/sqrt(num_somata)*(bounds(2,2)-bounds(1,2)) + bounds(1,2));
% %         l_info{l}.soma_loc = [zeros(num_somata,1), x(:), y(:)];
%         
%         l_info{l}.contact_loc = rand(1000000, 3);
%         l_info{l}.contact_loc(:,2) = l_info{l}.contact_loc(:,2) * ...
%             (valid_region(2,1) - valid_region(1,1)) + valid_region(1,1);
%         l_info{l}.contact_loc(:,3) = l_info{l}.contact_loc(:,3) * ...
%             (valid_region(2,2) - valid_region(1,2)) + valid_region(1,2);
% 
%         l_info{l}.cell_ids = ceil(rand(1000000,2) * max(l_info{l}.cell_ids(:)));
% 
%         valid_contacts = l_info{l}.cell_ids(:,1) ~= l_info{l}.cell_ids(:,2);
% 
%         soma_locs = [l_info{l}.soma_loc(l_info{l}.cell_ids(:,1),2:3), ...
%             l_info{l}.soma_loc(l_info{l}.cell_ids(:,2),2:3)];
% 
%         rel_locs = [l_info{l}.contact_loc(:,2:3), l_info{l}.contact_loc(:,2:3)] - soma_locs;
%         dists = [sqrt(sum(rel_locs(:,1:2).^2,2)), sqrt(sum(rel_locs(:,3:4).^2,2))];
%         
%         rel_locs(:,1) = rel_locs(:,1)./dists(:,1);
%         rel_locs(:,2) = rel_locs(:,2)./dists(:,1);
%         rel_locs(:,3) = rel_locs(:,3)./dists(:,2);
%         rel_locs(:,4) = rel_locs(:,4)./dists(:,2);
%         angles = acos(rel_locs(:,1).*rel_locs(:,3) + rel_locs(:,2).*rel_locs(:,4));        
% 
%         valid_contacts = valid_contacts & dists(:,1) < dist_thresh & dists(:,2) < dist_thresh;
%         angles = angles(valid_contacts);
% 
% 
%         layer_plot_data{l} = hist(angles, angle_bins);
%         
% %         figure; scatter(l_info{l}.soma_loc(:,2), l_info{l}.soma_loc(:,3), 'fill');
% %         scatter(soma_locs(1:100:end,2), soma_locs(1:100:end,3), 'fill');
% %         hold all; scatter(l_info{l}.contact_loc(1:100:end,2), l_info{l}.contact_loc(1:100:end,3), 'fill');
% %         scatter(soma_locs(valid_contacts,2), soma_locs(valid_contacts,3), 'fill');
%     end
%         
%          
%     figure; hold all
%     plot(angle_bins, layer_plot_data{1}/sum(layer_plot_data{1}), 'lineWidth', 2);
%     plot(angle_bins, layer_plot_data{2}/sum(layer_plot_data{2}), 'lineWidth', 2);
%     
% 
%     cns =  {[70014, 70158, 70016, 70154], [70164, 70176, 70230, 70205]};
% %     cns =  {[70130, 70080, 70147, 70116], [70197, 70169, 70244, 70191]};
%     legend_lbl = cell(2,1);
%     
%     for l = 1:2
%         
%         figure; hold all
%         for n = 1:length(cns{l})
%             legend_lbl{l}{n} = num2str(cns{l}(n));
%             l_info{l} = contact_info(cns{l}(n), C.type.(l_names{l}));
% %         l_info{l} = contact_info(70014, C.type.(l_names{l}));
% %         l_info{l} = contact_info(70158, C.type.(l_names{l}));
% %         valid_somata = ...
% %             l_info{l}.soma_loc(:,2) > valid_region(1,1) & ...
% %             l_info{l}.soma_loc(:,3) > valid_region(1,2) & ...
% %             l_info{l}.soma_loc(:,2) < valid_region(2,1) & ...
% %             l_info{l}.soma_loc(:,3) < valid_region(2,2);
%         
% 
%             valid_contacts = l_info{l}.cell_ids(:,1) ~= l_info{l}.cell_ids(:,2);
%             valid_contacts = valid_contacts & l_info{l}.contact_size > size_threshold;
%             valid_contacts = valid_contacts & ~(...
%                 l_info{l}.contact_loc(:,2) > valid_region(1,1) & ...
%                 l_info{l}.contact_loc(:,3) > valid_region(1,2) & ...
%                 l_info{l}.contact_loc(:,2) < valid_region(2,1) & ...
%                 l_info{l}.contact_loc(:,3) < valid_region(2,2));
% %         valid_contacts = valid_contacts & valid_somata(l_info{l}.cell_ids(:,1)) & ...
% %             valid_somata(l_info{l}.cell_ids(:,2));
%                                                                
%             for k = 1:length(list_names);
%                 l_info{l}.(list_names{k}) = l_info{l}.(list_names{k})(valid_contacts,:);
%             end
% 
%             soma_locs = [l_info{l}.soma_loc(l_info{l}.cell_ids(:,1),2:3), ...
%                 l_info{l}.soma_loc(l_info{l}.cell_ids(:,2),2:3)];
% 
%             rel_locs = [l_info{l}.contact_loc(:,2:3), l_info{l}.contact_loc(:,2:3)] - soma_locs;
%             dists = [sqrt(sum(rel_locs(:,1:2).^2,2)), sqrt(sum(rel_locs(:,3:4).^2,2))];
% 
%             rel_locs(:,1) = rel_locs(:,1)./dists(:,1);
%             rel_locs(:,2) = rel_locs(:,2)./dists(:,1);
%             rel_locs(:,3) = rel_locs(:,3)./dists(:,2);
%             rel_locs(:,4) = rel_locs(:,4)./dists(:,2);
%             angles = acos(rel_locs(:,1).*rel_locs(:,3) + rel_locs(:,2).*rel_locs(:,4));   
% 
%             layer_plot_data{l} = hist(angles, angle_bins);
% 
%             plot(angle_bins, layer_plot_data{l}/sum(layer_plot_data{l}), 'lineWidth', 2);
%         end
%         legend(legend_lbl{l})
%     end
end
        