function sac_to_sac_connections

C = get_constants;
% dsgc_nums = C.type.ganglion;
% dsgc_nums = C.type.oodsgc;
% dsgc_nums = C.type.on_dsgc;
% wide_oo = [20221 20187 17069];
% dsgc_nums = wide_oo;

% dsgc_nums = [[20181 17140 20208 20178 17097 17114 17084 20140 20129 30003 20071 30002 20019 20011 20005 26057], ...
%     [17216 10005 10013 26062 26118], ...
%     [50004 17013 17092 20051 20234 20082]];
% dsgc_nums = 20082;
gc_types = cell_info_typedef_gc;

sac_nums{1} = C.type.sure_off_sac;
sac_nums{2} = C.type.on_sac;

max_cont_size = 500;

num_bins = 24;
vericose_only = true;

nhood_size = 1000;
p = cell(2,1);


fns = fieldnames(gc_types);
for t = 32:length(fns);
    dsgc_nums = gc_types.(fns{t});



for dn = 1:length(dsgc_nums);
    try
    d = dsgc_nums(dn);
    
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
                
                for k = 1:length(cont_angles)
                    kn = ceil(180*(cont_angles(k)/pi + 1));
                    angle_count(kn, l) = angle_count(kn, l) + 1;
                    angle_num(kn, l) = angle_num(kn, l) + min(s_conts(2,k), max_cont_size);
                    
                end
            end
        end
       
    
    end
    
    
    %rebin
    
    count_plot = rebin_data(angle_count, num_bins);        
    cont_size_plot = rebin_data(angle_num, num_bins)./count_plot;        
    area_normed_plot = rebin_data(angle_num, num_bins)./rebin_data(angle_denom, num_bins);         
    count_normed_plot = rebin_data(angle_count, num_bins)./rebin_data(angle_denom, num_bins);
    denom_plot = rebin_data(angle_denom, num_bins);
    
    plot_theta = (([(1:num_bins)'; 1]-.5) * 2 * pi / num_bins) * ones(1,size(angle_num,2)) + pi;
    
    
    
    figure;
    
    subplot(2,2,1);       
    polar(plot_theta, area_normed_plot);
    title('normalized contact area');
    
    subplot(2,2,2);
    polar(plot_theta, count_normed_plot);
    title('normalized contact count');
    
    subplot(2,2,3);
    polar(plot_theta, cont_size_plot);
    title('mean contact size');
    
    subplot(2,2,4);   
    polar(plot_theta, count_plot);
    title('contact count');
    
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
    

    set(gcf, 'Position', [0 0 800 800]);
    
    try
        mkdir(['~/data/stratification/images/sac2gc/' fns{t}]);
    catch
    end
    
    saveas(gcf, ['~/data/stratification/images/sac2gc/' fns{t} '/sac2' num2str(d) '_' cell_type(d) '.png']);
%     saveas(gcf, ['~/data/stratification/images/sac2' num2str(d) '_' cell_type(d) '.png']);
    close all
    
     
    catch ME
        disp(ME.message)
    end
end

end
    
end
            
            
        

function new_data = rebin_data(data, num_bins)
    
    for k = 1:size(data,2)
        temp_data = reshape(data(:,k), [size(data,1)/num_bins, num_bins]);
        temp_data = sum(temp_data,1);
        temp_data = [temp_data'; temp_data(1)];
        new_data(:,k) = temp_data;
    end
    
    
end
        
        
        
    
    
    