function sac_to_sac_connections

close all

C = get_constants;
% dsgc_nums = C.type.ganglion;
% dsgc_nums = C.type.oodsgc;
dsgc_nums = {[90001 17080 20213 25005 20220], [90002 20125], [20239 20254 20245 20179 20210], [20233 17161 20137 20096]};

distance_hist = cell(length(dsgc_nums),1);

% dsgc_nums = C.type.on_dsgc;
% wide_oo = [20221 20187 17069];
% dsgc_nums = wide_oo;
sac_nums{1} = C.type.sure_off_sac;
sac_nums{2} = C.type.on_sac;

max_cont_size = 500;

num_bins = 24;
vericose_only = true;

nhood_size = 1000;
p = cell(2,1);

plot_data = cell(2, length(dsgc_nums));

for dnt = 1:length(dsgc_nums);
    
    distance_hist{dnt} = zeros(250, length(dsgc_nums{dnt}), 2, 4);
    
%     fh(1) = figure; hold all;    
%     polar([], []);
%     ax_h(1) = gca;
%     
%     fh(2) = figure; hold all;
%     polar([], []);
%     ax_h(2) = gca;
    
for dn = 1:length(dsgc_nums{dnt});
    try
    d = dsgc_nums{dnt}(dn);
    
    c_d = cell_data(d);
    gc_mid = c_d.get_midpoint(true);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
    p{1} = d_p(p_depth > 10 & p_depth < 50, 2:3);
    p{2} = d_p(p_depth > 50 & p_depth < 80, 2:3);
    
%     d_p = [];
    
    angle_denom = zeros(360,2);
    angle_num = zeros(360,2);
%     angle_size_total = zeros(360,2);
    angle_count = zeros(360,2);
    
    hull = cell(2,1);
    
    for l = 1:2
%         if l == 1
%             layer_name = 'off';
%         else
%             layer_name = 'on';
%         end
        
        
        p_min = min(p{l});
%         p_max = max(p{l});
        

        hull_inds = convhull(p{l}(:,1),p{l}(:,2));
        hull{l} = [p{l}(hull_inds,1) p{l}(hull_inds,2)];

        for k = 1:2
            p{l}(:,k) = ceil((p{l}(:,k) - p_min(k) + nhood_size/2)/nhood_size);
        end
        
        
        
        
        locality_mat = false(max(p{l}(:,1)), max(p{l}(:,2)));
        
        locality_mat(sub2ind(size(locality_mat), p{l}(:,1), p{l}(:,2))) = true;
        
        for s = 1:length(sac_nums{l})
            
            c_d = cell_data(sac_nums{l}(s));
            
            s_mid = c_d.get_midpoint;
            
            sp = c_d.get_surface;
            sp_depth = C.f(sp(:,1));
            
            
            sp = sp(sp_depth > 15 & sp_depth < 85, 2:3);
            angles = atan2(sp(:,2)-s_mid(3), sp(:,1)-s_mid(2));
            
            for k = 1:2
                sp(:,k) = ceil((sp(:,k) - p_min(k) + nhood_size/2)/nhood_size);
            end
            
            is_valid = sp(:,1) > 0 & sp(:,2) > 0 & sp(:,1) < size(locality_mat,1) & sp(:,2) < size(locality_mat,2);
            
            sp = sp(is_valid,:);
            angles = angles(is_valid);
            
            is_near_dsgc = locality_mat(...
                sub2ind(size(locality_mat), ...
                sp(:,1), sp(:,2)));
            
            angles = angles(is_near_dsgc);
            angles = ceil((angles / pi + 1) * 180);
            
            angle_denom(:,l) = angle_denom(:,l) + hist(angles, 1:360)';
            
            
            if vericose_only
                s_conts = detect_vericose_contacts(sac_nums{l}(s), 500, 200, 50000);
            else
                s_conts = double(c_d.contacts);
            end
            
            
            if ~isempty(s_conts)
                
                s_conts = s_conts(:, s_conts(1,:) == d);
            
                cont_angles = atan2(s_conts(5,:)-s_mid(3), s_conts(4,:)-s_mid(2));
                cont_dist = ceil(sqrt( (s_conts(4,:) - s_mid(2)).^2 + (s_conts(5,:) - s_mid(3)).^2 ) / 1000);                
                angle_bin = ( cont_angles < pi/4 & cont_angles >= -pi/4 ) + ...
                    2 * ( cont_angles < 3 * pi/4 & cont_angles >= pi/4 ) + ...
                    3 * ( cont_angles >= 3 * pi/4 | cont_angles < -3 * pi/4 ) + ...
                    4 * ( cont_angles >= -3 * pi/4 & cont_angles < -pi/4 );
                
                for k = 1:length(cont_angles)
                    
                    kn = ceil(180*(cont_angles(k)/pi + 1));
                    angle_count(kn, l) = angle_count(kn, l) + 1;
                    angle_num(kn, l) = angle_num(kn, l) + min(s_conts(2,k), max_cont_size);                                        
                    
                    distance_hist{dnt}(cont_dist(k), dn, l, angle_bin(k)) = ...
                        distance_hist{dnt}(cont_dist(k), dn, l, angle_bin(k)) + 1;
                end
            end
        end
       
    
    end
    
    
    %rebin
    
    area_normed_plot = rebin_data(angle_num, num_bins)./rebin_data(angle_denom, num_bins);         
    
    plot_data{1, dnt} = [plot_data{1, dnt}, area_normed_plot(:,1)];
    plot_data{2, dnt} = [plot_data{2, dnt}, area_normed_plot(:,2)];
    
    
    
    
    
%     figure;
    
%     subplot(2,2,1);       
        
%     polar(ax_h(1), plot_theta(:,1), );
%     title(num2str(dnt));
%     polar(ax_h(2), plot_theta(:,1), area_normed_plot(:,2));
%     title(num2str(dnt));
    
    
%     subplot(2,2,2);
%     polar(plot_theta, count_normed_plot);
%     title('normalized contact count');
%     
%     subplot(2,2,3);
%     polar(plot_theta, cont_size_plot);
%     title('mean contact size');
%     
%     subplot(2,2,4);   
%     polar(plot_theta, count_plot);
%     title('contact count');
    
    % hold all;
%     plot(hull{1}(:,1), hull{1}(:,2), 'lineWidth', 2);
%     plot(hull{2}(:,1), hull{2}(:,2), 'lineWidth', 2);
%     
%     plot([C.min_xy(2), C.max_xy(2), C.max_xy(2), C.min_xy(2), C.min_xy(2),]', ...
%         [C.min_xy(1), C.min_xy(1), C.max_xy(1), C.max_xy(1), C.min_xy(1)]', ...
%         'Color', [0 0 0], 'lineWidth', 2);
%     
%     scatter(gc_mid(2), gc_mid(3), '*', 'markerEdgeColor', [0 0 0]);
%     
%     title([num2str(d) ' hull']);
    

%     set(gcf, 'Position', [0 0 800 800]);
%     saveas(gcf, ['~/data/stratification/images/sac2' num2str(d) '_' cell_type(d) '.png']);
%     close all
    
     
    catch ME
        disp(ME.message)
    end
end

% 
%     plot_theta = (([(1:num_bins)'; 1]-.5) * 2 * pi / num_bins) * ones(1,size(plot_data{1, dnt},2)) + pi;
%     figure; polar(plot_theta, plot_data{1, dnt}); 
%     hs = get(gca, 'Children');
%     for k = hs'; set(k, 'lineWidth', 2); end
%     save_fig([C.image_dir, 'oodsgc_off_conn_' num2str(dnt) '.eps'], gcf);
%     
%     figure; polar(plot_theta, plot_data{2, dnt});
%     hs = get(gca, 'Children');
%     for k = hs'; set(k, 'lineWidth', 2); end
%     save_fig([C.image_dir, 'oodsgc_on_conn_' num2str(dnt) '.eps'], gcf);
%     
    
end
    

theta = (([(1:num_bins)'; 1]-.5) * 2 * pi / num_bins) * ones(1,size(plot_data{1, dnt},2)) + pi;
% 
% M = zeros(4,4,2);
% for n = 1:4
%     for d = 1:4; 
%         for l = 1:2;
%                 a = squeeze(distance_hist{n}(:,:,l,d));
%                 for x = 1:size(a,2); 
%                     a(:,x) = a(:,x)/sum(a(:,x)); 
%                 end
%                 disp(sum(a));
%                 a = a'*(1:size(a,1))';
%                 M(n,d,l) = mean(a);
%         end
%     end 
% end

save('~/data/stratification/model_dump.mat', 'distance_hist', 'plot_data');

dist_hist = zeros(size(distance_hist{1},1), 2);

for l = 1:2
    for n = 1:4
        t = squeeze(distance_hist{n}(:,:,l,:));
        dist_hist(:,l) = dist_hist(:,l) + sum(t(:,:), 2);
    end
    
    K = gausswin(11);
    K = K/sum(K);
    dist_hist(:,l) = dist_hist(:,l) / sum(dist_hist(:,l));
    dist_hist(:,l) = conv(dist_hist(:,l), K, 'same');
end

for l = 1:2;    
    
    rot_plot_data{l} = plot_data{l, 1}(1:end-1,:);
    rot_plot_data{l} = [rot_plot_data{l}, ...
        plot_data{l, 2}([3*num_bins/4+1:end-1, 1:3*num_bins/4], :)];
    rot_plot_data{l} = [rot_plot_data{l}, ...
        plot_data{l, 3}([2*num_bins/4+1:end-1, 1:2*num_bins/4], :)];
    rot_plot_data{l} = [rot_plot_data{l}, ...
        plot_data{l, 4}([num_bins/4+1:end-1, 1:num_bins/4], :)];
end
        
angle_distribution = [mean(rot_plot_data{1}, 2), mean(rot_plot_data{2}, 2)];


end
            
            
        

function new_data = rebin_data(data, num_bins)
    
    for k = 1:size(data,2)
        temp_data = reshape(data(:,k), [size(data,1)/num_bins, num_bins]);
        temp_data = sum(temp_data,1);
        temp_data = [temp_data'; temp_data(1)];
        new_data(:,k) = temp_data;
    end
    
    
end
        
        
        
    
    
    