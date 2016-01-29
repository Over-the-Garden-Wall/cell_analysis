function sac_to_sac_connections

C = get_constants;
dsgc_nums = [17152 17053 20180 17069 20251 20246 20250 20241];

sac_nums = C.type.on_sac;

max_cont_size = 500;

num_bins = 24;
    


nhood_size = 1000;
p = cell(2,1);

for dn = 1:length(dsgc_nums);
    
    d = dsgc_nums(dn);
    
    c_d = cell_data(d);
    gc_mid = c_d.get_midpoint(true);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
%     p{1} = d_p(p_depth > 10 & p_depth < 50, 2:3);
    p = d_p(p_depth > 50 & p_depth < 80, 2:3);
    
%     d_p = [];
    
    angle_denom = zeros(360,1);
    angle_num = zeros(360,1);
    angle_count = zeros(360,1);
    
%     hull = cell(2,1);
    
%     for l = 1:2
        
        p_min = min(p);
%         p_max = max(p);
        

        hull_inds = convhull(p(:,1),p(:,2));
        hull = [p(hull_inds,1) p(hull_inds,2)];

        for k = 1:2
            p(:,k) = ceil((p(:,k) - p_min(k) + nhood_size/2)/nhood_size);
        end
        
        
        
        
        locality_mat = false(max(p(:,1)), max(p(:,2)));
        
        locality_mat(sub2ind(size(locality_mat), p(:,1), p(:,2))) = true;
        
        for s = 1:length(sac_nums)
            
            c_d = cell_data(sac_nums(s));
            
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
            
            angle_denom = angle_denom + hist(angles, 1:360)';
            
            
            
            s_conts = double(c_d.contacts);
            s_conts = s_conts(:, s_conts(1,:) == d);
            
            
            cont_angles = atan2(s_conts(5,:)-s_mid(3), s_conts(4,:)-s_mid(2));
            for k = 1:length(cont_angles)
                kn = ceil(180*(cont_angles(k)/pi + 1));
                angle_count(kn) = angle_count(kn) + 1;
                angle_num(kn) = angle_num(kn) + min(s_conts(2,k), max_cont_size);
            end
        end
        
    
    
    
    %rebin
    plot_data = angle_num./angle_denom;
    plot_data = reshape(plot_data, [num_bins, 360/num_bins]);
    plot_data = squeeze(sum(plot_data,1));
    plot_data = [plot_data, plot_data(1,:)];
    
    plot_theta = ((1:360)'/180-1)*pi;
    plot_theta = reshape(plot_theta, [num_bins, 360/num_bins]);
    plot_theta = squeeze(mean(plot_theta,1));
    plot_theta = [plot_theta, plot_theta(1,:)];
    
    figure;
    
    subplot(2,2,1);
    
    
    polar(plot_theta, plot_data);
    title(num2str(d));
    
    
    num_bins = 24;
    plot_data = angle_denom;
    plot_data = reshape(plot_data, [num_bins, 360/num_bins]);
    plot_data = squeeze(sum(plot_data,1));    
    plot_data = [plot_data, plot_data(1,:)];
    
    
    subplot(2,2,2);
    polar(plot_theta, plot_data);
    title([num2str(d) ' denom']);
    
    subplot(2,2,3); hold all
    
    plot(hull(:,1), hull(:,2), 'lineWidth', 2);
    
    plot([C.min_xy(2), C.max_xy(2), C.max_xy(2), C.min_xy(2), C.min_xy(2),]', ...
        [C.min_xy(1), C.min_xy(1), C.max_xy(1), C.max_xy(1), C.min_xy(1)]', ...
        'Color', [0 0 0], 'lineWidth', 2);
    
    scatter(gc_mid(2), gc_mid(3), '*', 'markerEdgeColor', [0 0 0]);
    
    title([num2str(d) ' hull']);
    
    set(gcf, 'Position', [0 0 800 800]);
    saveas(gcf, ['~/data/stratification/images/sac2ONdsgc' num2str(d) '.png']);
    close all
end
    
end
            
            
        
        
        
        
    
    
    