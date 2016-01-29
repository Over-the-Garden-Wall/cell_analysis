function dsgc_rf2(gc, path_length, min_dist)

    C = get_constants;

    bip_radius = 10000;
    pix_resolution = 500;
    
    sacs{1} = C.type.sure_off_sac;
    sacs{2} = C.type.on_sac;
    
%     type_switch = [0, 50000, 90000; 0, 30000, 90000]; 
    type_switch = [45000, 50000, 55000; 25000, 30000, 35000];
    
    c_d = cell_data(gc);
    gc_mid = mean(c_d.hull_2d);
    
%     target_region = [];
    target_region = [gc_mid - 25000; gc_mid + 25000];
%     target_region = [gc_mid + 35000; gc_mid + 65000];
%     target_region(:,1) = target_region(:,1) + [-1; 1]*150000;
    
    layer_depths = [0 40; 50 80];
    
    for l = 1:2
        type_com = zeros(2,3);
        type_points = cell(2,1);
    
        all_vconns = [];
        
        for n = 1:length(sacs{l})
            conns = detect_vericose_contacts(sacs{l}(n), 500, 3, 1, 60000);
            
            is_me = conns(1,:) == gc;
            if any(is_me)
                sub_conns = conns(:,is_me);
                sub_conns(:,sub_conns(2,:) < 100) = [];
                all_vconns = [all_vconns, [sacs{l}(n) * ones(1,size(sub_conns,2)); sub_conns]];
            end
            
            if ~isempty(target_region)
            is_valid = all_vconns(5,:) >= target_region(1,1) & ...
                all_vconns(5,:) <= target_region(2,1) & ...
                all_vconns(6,:) >= target_region(1,2) & ...
                all_vconns(6,:) <= target_region(2,2);
            all_vconns = all_vconns(:,is_valid);
            end
            
        end
        
        sac_paths = cell(size(all_vconns,2),1);
        
        
        for n = 1:size(all_vconns,2)
            
            sac_paths{n} = zeros(path_length, 3);
            
            s = load([C.skele_dir 's' num2str(all_vconns(1,n)) '.mat']);
            
            sc_d = cell_data(all_vconns(1,n));
            
            d = sqrt(sum((s.nodes - ones(size(s.nodes,1),1) * (all_vconns(4:6,n))').^2,2));
            [dummy, minind] = min(d);
            
            sac_paths{n}(end,:) = s.nodes(minind,:);
                
            sac_soma = sc_d.get_midpoint;
            
            for k = path_length-1:-1:1
                next_node = sac_paths{n}(k+1,:);
            
                dist_to_soma = sqrt(sum((sac_soma - next_node).^2));
                
                step_size = (dist_to_soma - min_dist) / k;
                step_axis = (sac_soma(2:3) - next_node(2:3));
                step_axis = step_axis / sqrt(sum(step_axis.^2));
                
                next_node(2:3) = next_node(2:3) + step_size * step_axis;
                
                d = sqrt(sum((s.nodes - ones(size(s.nodes,1),1) * next_node).^2,2));
                [dummy, minind] = min(d);
                sac_paths{n}(k,:) = s.nodes(minind,:);
                                            
            end
            
            d = sqrt(sum((sac_paths{n} - ones(size(sac_paths{n},1),1) * sac_soma).^2,2));
            
            for k = 1:2
                type_points{k} = [type_points{k}; sac_paths{n}(d > type_switch(l,k) & d < type_switch(l,k+1),:)];                
            end
            
        end
        
        for k = 1:2
            for p = 1:size(type_points{k})
                type_com(k,1:2) = type_com(k,1:2) + type_points{k}(p, 2:3);
                type_com(k,3) = type_com(k,3) + 1;
            end
            type_com(k,1:2) = type_com(k,1:2) / type_com(k,3);
        end
                
        
        min_points = min([type_points{1}; type_points{2}]);
        max_points = max([type_points{1}; type_points{2}]);
        
        pix_map = zeros([ceil((max_points(2:3) - min_points(2:3) + 1)/pix_resolution), 2]);
        
        for k = 1:2
            pix_points = [type_points{k}(:,2) - min_points(2), type_points{k}(:,3) - min_points(3)];
            pix_points = ceil((pix_points+1)/pix_resolution);
            
            for n = 1:size(pix_points,1)
                pix_map(pix_points(n,1), pix_points(n,2), k) = pix_map(pix_points(n,1), pix_points(n,2), k) + 1;
            end
        end
        
        k_bound = ceil(bip_radius/pix_resolution);
        [x, y] = meshgrid(-k_bound:k_bound, -k_bound:k_bound);
        r = sqrt(x.^2 + y.^2);
        K = double(r*pix_resolution <= bip_radius);               
        
        type_contact_image = zeros(size(pix_map) + [size(K)-1, 0]);
        for k = 1:2
            type_contact_image(:,:,k) = conv2(pix_map(:,:,k), K, 'full');
        end
        type_contact_image = permute(type_contact_image, [3 1 2]);
        type_contact_image = log(type_contact_image+1);
        type_contact_image(1,:) = type_contact_image(1,:) / max(type_contact_image(1,:));
        type_contact_image(2,:) = type_contact_image(2,:) / max(type_contact_image(2,:));
        
        im_colors = [-.5 -.5 0; 0 -.5 -.5]' * type_contact_image(:,:);
        im_colors = im_colors + 1;
        
        [X, Y] = meshgrid(1:size(type_contact_image,3), 1:size(type_contact_image,2));
        Y = Y*pix_resolution + min_points(2) - bip_radius + pix_resolution/2;
        X = X*pix_resolution + min_points(3) - bip_radius + pix_resolution/2;
        
        figure; hold all;
        for n = 1:size(im_colors, 2)
            if ~all(im_colors(:,n)==1)
                rectangle('Position', [Y(n), X(n), pix_resolution, pix_resolution], ...
                    'edgeColor', im_colors(:,n)', 'faceColor', im_colors(:,n)');
            end
        end
%         plot_cells(gc, 1, .02, [.5 .5 .5]);
        
        p = c_d.get_surface;
        d = C.f(p(:,1));
        p = p(d >= layer_depths(l,1) & d <= layer_depths(l,2), 2:3);
        mean_p = mean(p);
        
        is_hull = convhull(p(:,1), p(:,2));
        h = cell(1,2);
        [h{:}] = poly2cw(p(is_hull,1), p(is_hull,2));
        
        plot(h{1}, h{2}, 'k', 'lineWidth', 2);
        
        scatter(type_com(1,1), type_com(1,2), 50, 'markerFaceColor', [0 0 1]);
        scatter(type_com(2,1), type_com(2,2), 50, 'markerFaceColor', [1 0 0]);
        
        if ~isempty(target_region)
            rectangle('Position', [min(target_region(:,1)), min(target_region(:,2)), ...
                abs(target_region(2,1) - target_region(1,1)), abs(target_region(2,2) - target_region(1,2))], ...
                        'edgeColor', [0 0 0]);
        else
            scatter(mean_p(1), mean_p(2), 50, 'markerFaceColor', [0 0 0]);        
        end
                
                
%         h = cell(2,2);
%         crc = {cos((0:360)/180*pi)'*bip_radius, sin((0:360)/180*pi)'*bip_radius};
%         [crc{:}] = poly2cw(crc{1}, crc{2});
% %         for k = 1:2
% %             for p = 1:size(type_points{k})
% %                 [h{k,:}] = polybool('union', h{k,1}, h{k,2}, crc{1} + type_points{k}(p,2), crc{2} + type_points{k}(p,3));
% %             end
% %         end
%         for k = 1:2
%             while ~isempty(type_points{k})
%                 [h{k,:}] = polybool('union', h{k,1}, h{k,2}, crc{1} + type_points{k}(1,2), crc{2} + type_points{k}(1,3));
%               d = sqrt((type_points{k}(:,2) - type_points{k}(1,2)).^2 + (type_points{k}(:,3) - type_points{k}(1,3)).^2);
%               type_points{k}(d < bip_radius/2,:) = [];
%             end
%         end
%         figure; hold all;
%         plot(h{1,1}, h{1,2});
%         plot(h{2,1}, h{2,2});
        
    end