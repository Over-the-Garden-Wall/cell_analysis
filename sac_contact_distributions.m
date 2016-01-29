function sac_contact_distributions

    C = get_constants;
    
    sac_type = 'on_sac';
    num_angle_bins = 36;
    angle_span = pi / 6;
    min_coverage = .8;
    sac_input_length = 90000;
    num_dp_to_visualize = 5;
    
    
    line_points = 100;
    
    if strcmp(sac_type, 'on_sac')
        types{1} = {'BC7'};
        types{2} = {'BC5t', 'BC5o', 'BC5i'};
    else
        types{1} = {'BC1', 'BC2'};
        types{2} = {'BC3a', 'BC3b'};
    end
    
    sacs = C.type.(sac_type);
    num_sacs = length(sacs);
    num_tsets = length(types);
    
    ts_hull = cell(num_tsets,2);
    for ts = 1:num_tsets
        for t = 1:length(types{ts})
            for c = C.type.(types{ts}{t})
                c_d = cell_data(c);
                [h{1}, h{2}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
                [ts_hull{ts,:}] = polybool('union', h{1}, h{2}, ts_hull{ts,1}, ts_hull{ts,2});
            end
        end
    end
    intersect_hull = cell(1,2);
    [intersect_hull{:}] = polybool('intersection', ts_hull{1,1}, ts_hull{1,2}, ts_hull{2,1}, ts_hull{2,2});
    
    angles = (.5:num_angle_bins)/num_angle_bins * 2 * pi;
    angle_bin_width = angles(2)-angles(1);
    
    contact_distance = cell(1, num_tsets);
    contact_angle = cell(1, num_tsets);
    contact_weight = cell(1, num_tsets);
            
    for sac_n = 1:num_sacs
        sac = sacs(sac_n);
        c_d = cell_data(sac);
        sac_mid = c_d.get_midpoint(true);
        sac_conts = double(c_d.contacts);
        
        
        for ts = 1:num_tsets
            is_my_type = false(99999,1);
            
            for t = 1:length(types{ts})
                is_my_type(C.type.(types{ts}{t})) = true;
            end
            my_conts = sac_conts(:,is_my_type(sac_conts(1,:)));
            
            contact_distance{ts} = sqrt((my_conts(4,:) - sac_mid(2)).^2 + (my_conts(5,:) - sac_mid(3)).^2);
            contact_angle{ts} = atan2((my_conts(5,:) - sac_mid(3)), (my_conts(4,:) - sac_mid(2)));
            contact_angle{ts}(contact_angle{ts}<0) = contact_angle{ts}(contact_angle{ts}<0) + 2*pi;            
            contact_weight{ts} = my_conts(2,:);
%             cont_mean(sac_n,ts) = sum(cont_weight .* dist_from_soma);
%             cont_msd = (dist_from_soma - cont_mean(sac_n,ts)).^2;
%             cont_sd(sac_n,ts) = sum(cont_weight .* cont_msd);
        end

        cont_mean = nan(num_angle_bins, num_tsets);
        cont_sd = nan(num_angle_bins, num_tsets);
        
        
        for theta_n = 1:num_angle_bins
            theta = angles(theta_n);
            
            line = [sac_mid(2) + sin(theta)*(1:line_points)'/line_points * sac_input_length, ...
                sac_mid(3) + cos(theta)*(1:line_points)'/line_points * sac_input_length];
            
            proportion_in_hull = sum(inpolygon(line(:,1), line(:,2), intersect_hull{1}, intersect_hull{2}))/line_points;
            
            
            
            if proportion_in_hull > min_coverage
                
                for ts = 1:num_tsets
                    adjusted_angles = contact_angle{ts};
                    if theta < pi
                        adjusted_angles(adjusted_angles>pi*1.5) = adjusted_angles(adjusted_angles>pi*1.5) - 2*pi;
                    else
                        adjusted_angles(adjusted_angles<pi*.5) = adjusted_angles(adjusted_angles<pi*.5) + 2*pi;
                    end
                    valid_contacts = adjusted_angles > theta - angle_span/2 & adjusted_angles < theta + angle_span/2;
                    c_weight = contact_weight{ts}(valid_contacts);
                    c_weight = c_weight/sum(c_weight);
                    c_dist = contact_distance{ts}(valid_contacts);
                    
                    
                    cont_mean(theta_n,ts) = sum(c_weight .* c_dist);
                    cont_msd = (c_dist - cont_mean(theta_n,ts)).^2;
                    cont_sd(theta_n,ts) = sum(c_weight .* cont_msd);                    
                    
                end
                
            end
        
        end
        
        %visualize!
        if sum(~isnan(cont_mean(:))) > num_dp_to_visualize

            figure; hold all;
            disjoint_polyfill(ts_hull{1,1}, ts_hull{1,2}, [1 0 0]);
            disjoint_polyfill(ts_hull{2,1}, ts_hull{2,2}, [0 0 1]);
            disjoint_polyfill(intersect_hull{1}, intersect_hull{2}, [1 0 1]);

            scatter(sac_mid(2), sac_mid(3), 30, '*', 'markerEdgeColor', [0 0 0]);

            colors = [1 1 0; .5 .5 .5; 0 0 0; 1 1 1];

            for ts = 1:num_tsets
                plot_data = [sac_mid(2) + sin(angles') .* cont_mean(:, ts), ...
                    sac_mid(3) + cos(angles') .* cont_mean(:, ts)];

                plot(plot_data(:,1), plot_data(:,2), ...
                    'lineWidth', 2,'Color', colors(ts,:), 'marker', '+', 'markerEdgeColor', colors(ts,:));
            end

            title(num2str(sac));
        
        end
    end
    
    
    
end
        
        
function disjoint_polyfill(X, Y, C)
    nan_locs = find(isnan(X));
    
    nan_locs = [0; nan_locs; length(X)+1];
    
    for n = 1:length(nan_locs)-1
        fill(X(nan_locs(n)+1:nan_locs(n+1)-1), Y(nan_locs(n)+1:nan_locs(n+1)-1), C);
    end
end
        
        