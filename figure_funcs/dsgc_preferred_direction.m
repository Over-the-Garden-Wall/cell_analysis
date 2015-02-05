function [mean_axis, pref_axes] = dsgc_preferred_direction(d)

C = get_constants;
sac_nums{1} = C.type.sure_off_sac;
sac_nums{2} = C.type.on_sac;

max_cont_size = 500;
    

nhood_size = 1000;
p = cell(2,1);
% 
% for dn = 1:length(dsgc_nums);
%     
%     d = dsgc_nums(dn);
    
    c_d = cell_data(d);
    gc_mid = c_d.get_midpoint(true);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
    p{1} = d_p(p_depth > 10 & p_depth < 50, 2:3);
    p{2} = d_p(p_depth > 50 & p_depth < 80, 2:3);
    
%     d_p = [];
    
    angle_denom = zeros(360,2);
    angle_num = zeros(360,2);
    angle_count = zeros(360,2);
    
    hull = cell(2,1);
    
    for l = 1:2
        
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
            
            
            
            s_conts = double(c_d.contacts);
            s_conts = s_conts(:, s_conts(1,:) == d);
            
            
            cont_angles = atan2(s_conts(5,:)-s_mid(3), s_conts(4,:)-s_mid(2));
            for k = 1:length(cont_angles)
                kn = ceil(180*(cont_angles(k)/pi + 1));
                angle_count(kn, l) = angle_count(kn, l) + 1;
                angle_num(kn, l) = angle_num(kn, l) + min(s_conts(2,k), max_cont_size);
            end
        end
        
    end
    
    norm_data = angle_num./angle_denom;
    norm_data = [norm_data(end-88:end,:); norm_data; norm_data(1:99,:)];
    kern = ones(180, 1)/180;
    
    smoothed_data = [conv(norm_data(:,1), kern, 'valid'), conv(norm_data(:,2), kern, 'valid')];
    [dummy, max_angles] = max(smoothed_data);
    
    [dummy, max_total] = max(sum(smoothed_data,2));
    
    max_angles = (max_angles/180-1)*pi;
    max_total = (max_total/180-1)*pi;
    
    pref_axes = [-cos(max_angles); -sin(max_angles)];
    mean_axis = [-cos(max_total); -sin(max_total)];
% end
    
end
            
            
        
        
        
        
    
    
    