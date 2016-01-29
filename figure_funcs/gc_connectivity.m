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

CELL_ID_MAX = 99999;
P_SPARSITY = 1000;

sac_nums{1} = C.type.sure_off_sac;
sac_nums{2} = C.type.on_sac;

max_cont_size = 500;

num_bins = 8;
vericose_only = true;

nhood_size = 1000;
p = cell(2,1);


gc_type_names = fieldnames(gc_types);
gc_nums = [];
for t = 1:length(gc_type_names);
    gc_nums = [gc_nums gc_types.(gc_type_names{t})];
end

gc_type_name = cell(length(gc_nums),1);
k = 0;
for t = 1:length(gc_type_names);
    for c = gc_types.(gc_type_names{t})
        k = k+1;
        gc_type_name{k} = gc_type_names{t};
    end
end


connecting_types = {'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5t', 'BC5i', 'BC5o', 'XBC', 'BC6', 'BC7', 'BC89', 'off_sac', 'on_sac'};



num_types = length(connecting_types);
num_gcs = length(gc_nums);

    type_hulls = cell(num_types, 2);
    h = cell(1,2);
    for k = 1:num_types
        for n = C.type.(connecting_types{k})
            c_d = cell_data(n);
            [h{1} h{2}] = poly2cw(c_d.hull_2d(:,1),c_d.hull_2d(:,2));
            [type_hulls{k,:}] = polybool('union', h{1}, h{2}, type_hulls{k,1}, type_hulls{k,2});
        end
    end
        
    num2type = zeros(CELL_ID_MAX, 1);
    for k = 1:num_types
        num2type(C.type.(connecting_types{k})) = k;
    end
    
    total_overlap = zeros([num_types num_gcs]); %pix1 in overlap, pix2 in overlap
    total_contact = zeros(size(total_overlap));
    
    sac_contact_angles = zeros(num_gcs, num_bins, 2);






for dn = 1:length(gc_nums);
    try
    d = gc_nums(dn);
    
    c_d = cell_data(d);
    gc_mid = c_d.get_midpoint(true);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
%     p{1} = d_p(p_depth > 10 & p_depth < 50, 2:3);
%     p{2} = d_p(p_depth > 50 & p_depth < 80, 2:3);
    
    d_p = d_p(p_depth > 10 & p_depth < 80, 2:3);
%     d_p = [];
    
    angle_denom = zeros(360,2);
    angle_num = zeros(360,2);
%     angle_size_total = zeros(360,2);
    angle_count = zeros(360,2);
    
%     hull = cell(2,1);
    
    p_min = min(d_p);
%         p_max = max(p{l});
        

%         hull_inds = convhull(p(:,1),p(:,2));
%         hull = [d_p(hull_inds,1) d_p(hull_inds,2)];

    
        for k = 1:2
            d_p(:,k) = ceil((d_p(:,k) - p_min(k) + nhood_size/2)/nhood_size);
        end
        
    for l = 1:2
%         if l == 1
%             layer_name = 'off';
%         else
%             layer_name = 'on';
%         end
        
        
        
        
        
        
        locality_mat = false(max(d_p(:,1)), max(d_p(:,2)));
        
        locality_mat(sub2ind(size(locality_mat), d_p(:,1), d_p(:,2))) = true;
        
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
    
    sac_contact_angles(dn, :, :) = area_normed_plot(1:end-1,:);
    
    figure;
    
    subplot(2,2,1);       
    polar(plot_theta, area_normed_plot);
    title('normalized contact area');
    
    subplot(2,2,2);   
    polar(plot_theta, count_plot);
    title('contact count');
    
%     subplot(2,2,3);
%     polar(plot_theta, cont_size_plot);
%     title('mean contact size');
    
    hand1 = subplot(2,2,3);
    
    for k = 1:num_types
            
                c_d = cell_data(gc_nums(dn));
                sp = c_d.get_surface;
                sp = sp(round(P_SPARSITY/2):P_SPARSITY:end,2:3);
                total_overlap(k,dn) = total_overlap(k,dn) + sum(inpolygon(sp(:,1), sp(:,2), type_hulls{k,1}, type_hulls{k,2}));
            
            
                conts = double(c_d.contacts);
                conts = conts(:,num2type(conts(1,:))==k);
                if ~isempty(conts)
                    total_contact(k,dn) = total_contact(k,dn) + sum(conts(2,:));
                end
            
        
        
    end
    
    
%     bar(3.5 + log(total_contact(:,dn)./M(:,dn)/P_SPARSITY)/log(10),1);
    bar(2+log(total_contact(:,dn)./total_overlap(:,dn)/P_SPARSITY * 100) / log(10));
    
%     set(gca, 'YScale', 'log');
            set(gca, 'YLim', [0 4], 'YTick', -0:4, 'YTickLabel', {'.01', '.1', '1', '10', '100'});                            
            set(gca, 'XTick', 1:num_types, 'XTickLabel', connecting_types);
            set(gca, 'fontSize', 20);
    
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
    
    hand2 = subplot(2,2,4);
    s = c_d.stratification;
    s = s(1-C.strat_x(1):end);
    s(end+1:101) = 0;
    s(102:end) = [];
    plot(s, 0:100, 'lineWidth', 2);
    set(gca, 'YDir', 'reverse')
    
    
    set(hand1, 'position', [.05 .1 .7 .4]);
    set(hand2, 'Position', [.8 .1, .15, .4]);

    set(gcf, 'Position', [0 0 800 1000]);
    
    try
        mkdir(['~/data/stratification/images/sac2gc/' gc_type_name{dn}]);
    catch
    end
    
    saveas(gcf, ['~/data/stratification/images/sac2gc/' gc_type_name{dn} '/sac2' num2str(d) '_' cell_type(d) '.png']);
%     saveas(gcf, ['~/data/stratification/images/sac2' num2str(d) '_' cell_type(d) '.png']);
    close all
    
    save('~/data/stratification/gc_contact_dump.mat', 'sac_contact_angles', 'total_overlap', 'total_contact');
    
     
    catch ME
        disp(ME.message)
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
        
        
        
    
    
    