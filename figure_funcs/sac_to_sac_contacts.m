function sac_to_sac_contacts

    C = get_constants;
    size_threshold = 0;
    num_angle_bins = 24;    
    
    angle_bins = (.5:num_angle_bins-.5) / num_angle_bins * pi;
    bin_step = angle_bins(2)-angle_bins(1);
    
    
    l_names = {'sure_off_sac', 'on_sac'};
    
    list_names = {'cell_ids', 'contact_loc', 'contact_size'};
    
    
    individual_cells = true;
%     cns =  [70014, 70158, 70016, 70154]; %off center
%     cns =  [70164, 70176, 70230, 70205]; %on center
cns = [70130, 70080, 70147, 70116]; %off l-edge
% cns = [70197, 70169, 70244, 70191]; %on l-edge
% cns = [70102 70134 70027 70089]; %off bl-corner
% cns = [70211 70215 70216 70188]; %on bl-corner

    if any(cns(1) == C.type.on_sac);
        l = 2;
    else
        l = 1;
    end
    for k = 1:length(cns)
        cns_lbl{k} = num2str(cns(k));
    end
%     cns =  {[70130, 70080, 70147, 70116], [70197, 70169, 70244, 70191]};

    
if ~individual_cells

    for l = 1:2
        l_info{l} = contact_info(C.type.(l_names{l}), C.type.(l_names{l}));
%         l_info{l} = contact_info(70014, C.type.(l_names{l}));

        num_sacs = size(l_info{l}.soma_loc,1);
        point_dist = zeros(250, num_sacs);
        for n = 1:num_sacs
            cn = l_info{l}.number(n);
            c_d = cell_data(cn);
            p = c_d.get_surface;
            d = C.f(p(:,1));
            p = p(d>20 & d<80,:);
            d = sqrt((p(:,2) - l_info{l}.soma_loc(n,2)).^2 + (p(:,3) - l_info{l}.soma_loc(n,3)).^2)/1000;
            point_dist(:,n) = hist(d, (1:250)-.5);% ./ (1:250);
            point_dist(:,n) = point_dist(:,n)/sum(point_dist(:,n));
        end

        valid_contacts = l_info{l}.number(l_info{l}.cell_ids(:,1)) ~= l_info{l}.number(l_info{l}.cell_ids(:,2));
        valid_contacts = valid_contacts & l_info{l}.contact_size > size_threshold;
                                                
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
        
        layer_num_data{l} = hist(angles, angle_bins);
        layer_den_data{l} = zeros(size(layer_num_data{l}));
        
        angle_weight = zeros(num_sacs,1);
        for n = 1:size(l_info{l}.contact_loc,1);
            cs = l_info{l}.cell_ids(n,:);
            
            rel_loc_to_somata = [l_info{l}.contact_loc(n,2) - l_info{l}.soma_loc(:,2), ...
                l_info{l}.contact_loc(n,3) - l_info{l}.soma_loc(:,3)];
            
            ds_to_somata = sqrt(sum(rel_loc_to_somata.^2,2));
            
            angles_to_somata = acos(rel_locs(n,1).*rel_loc_to_somata(:,1)./ds_to_somata + ...
                rel_locs(n,2).*rel_loc_to_somata(:,2)./ds_to_somata);
            
            for k = 1:num_sacs
                angle_weight(k) = point_dist(min(ceil(ds_to_somata(k)/1000),250), k);
            end
            angle_weight(cs(1)) = [];
            angle_weight = angle_weight/sum(angle_weight);
            angles_to_somata(cs(1)) = [];
            
            for k = 1:num_angle_bins
                is_this_angle = ...
                    angles_to_somata >= angle_bins(k)-bin_step/2 & ...
                    angles_to_somata < angle_bins(k)+bin_step/2;
                    
                layer_den_data{l}(k) = layer_den_data{l}(k) + sum(is_this_angle.*angle_weight);
                
            end
        end
        
        
    end
    
    figure; hold all
    plot(angle_bins, layer_num_data{1}./layer_den_data{1}, 'lineWidth', 2);
    plot(angle_bins, layer_num_data{2}./layer_den_data{2}, 'lineWidth', 2);
else
    figure; hold all;
    for cnk = 1:length(cns)
        l_info = contact_info(cns(cnk), C.type.(l_names{l}));
%         l_info{l} = contact_info(70014, C.type.(l_names{l}));

        num_sacs = size(l_info.soma_loc,1);
        point_dist = zeros(250, num_sacs);
        for n = 1:num_sacs
            cn = l_info.number(n);
            c_d = cell_data(cn);
            p = c_d.get_surface;
            d = C.f(p(:,1));
            p = p(d>20 & d<80,:);
            d = sqrt((p(:,2) - l_info.soma_loc(n,2)).^2 + (p(:,3) - l_info.soma_loc(n,3)).^2)/1000;
            point_dist(:,n) = hist(d, (1:250)-.5) ./ (1:250);
            point_dist(:,n) = point_dist(:,n)/sum(point_dist(:,n));
        end

        valid_contacts = l_info.number(l_info.cell_ids(:,1)) ~= l_info.number(l_info.cell_ids(:,2));
        valid_contacts = valid_contacts & l_info.contact_size > size_threshold;
                                                
        for k = 1:length(list_names);
            l_info.(list_names{k}) = l_info.(list_names{k})(valid_contacts,:);
        end
        
        soma_locs = [l_info.soma_loc(l_info.cell_ids(:,1),2:3), ...
            l_info.soma_loc(l_info.cell_ids(:,2),2:3)];
        
        rel_locs = [l_info.contact_loc(:,2:3), l_info.contact_loc(:,2:3)] - soma_locs;
        dists = [sqrt(sum(rel_locs(:,1:2).^2,2)), sqrt(sum(rel_locs(:,3:4).^2,2))];
        
        rel_locs(:,1) = rel_locs(:,1)./dists(:,1);
        rel_locs(:,2) = rel_locs(:,2)./dists(:,1);
        rel_locs(:,3) = rel_locs(:,3)./dists(:,2);
        rel_locs(:,4) = rel_locs(:,4)./dists(:,2);
        angles = acos(rel_locs(:,1).*rel_locs(:,3) + rel_locs(:,2).*rel_locs(:,4));   
        
        layer_num_data = hist(angles, angle_bins);
        layer_den_data = zeros(size(layer_num_data));
        
        angle_weight = zeros(num_sacs,1);
        for n = 1:size(l_info.contact_loc,1);
            cs = l_info.cell_ids(n,:);
            
            rel_loc_to_somata = [l_info.contact_loc(n,2) - l_info.soma_loc(:,2), ...
                l_info.contact_loc(n,3) - l_info.soma_loc(:,3)];
            
            ds_to_somata = sqrt(sum(rel_loc_to_somata.^2,2));
            
            angles_to_somata = real(acos(rel_locs(n,1).*rel_loc_to_somata(:,1)./ds_to_somata + ...
                rel_locs(n,2).*rel_loc_to_somata(:,2)./ds_to_somata));
            
            for k = 1:num_sacs
                angle_weight(k) = point_dist(min(ceil(ds_to_somata(k)/1000),250), k);
            end
            angle_weight(cs(1)) = [];
            angle_weight = angle_weight/sum(angle_weight);
            angles_to_somata(cs(1)) = [];
            
            for k = 1:num_angle_bins
                is_this_angle = ...
                    angles_to_somata >= angle_bins(k)-bin_step/2 & ...
                    angles_to_somata < angle_bins(k)+bin_step/2;
                    
                layer_den_data(k) = layer_den_data(k) + sum(is_this_angle.*angle_weight);
                
            end
        end
        
        plot(angle_bins, layer_num_data./layer_den_data, 'lineWidth', 2);
    end
%     
%     figure; hold all
%     
%     plot(angle_bins, layer_num_data{2}./layer_den_data{2}, 'lineWidth', 2);
    legend(cns_lbl);
    
end
end
        