function [observed_contacts contact_list] = contact_observedXdist(target_cells, contact_cells, max_dist, R)
    
    observed_contacts = zeros(max_dist,1);
%     bip_list = cell(length(target_cells),1);

    contact_list = [];

    for n = 1:length(target_cells)
        fn = ['./surface_points_trans/cell_' num2str(target_cells(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            
            [min_point, min_ind] = max(surface_points(:,1));
            
            soma_point = surface_points(min_ind,:);
        
            sub_Rs = pick_conns(R, target_cells(n), contact_cells);
            
            for k = 1:length(contact_cells)
                num_points = size(sub_Rs{k},2);
                if num_points > 0
                    dist = ceil(sqrt(sum((double(sub_Rs{k}(4:6,:)) - soma_point'*ones(1,num_points)).^2))/1000);                    
                    is_valid = dist <= max_dist;
                    if any(is_valid);                    
                        observed_contacts(dist(is_valid)) = ...
                            observed_contacts(dist(is_valid)) + double(sub_Rs{k}(3,is_valid))';                    

                    end
                    contact_list = [contact_list; dist' sub_Rs{k}(3,:)'];
%                     bip_list{k} = [bip_list{k} contact_cells(k)];
                end
            end
        end
    end
        
            
end