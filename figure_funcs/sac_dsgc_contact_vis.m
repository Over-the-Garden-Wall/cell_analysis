function sac_dsgc_contact_vis(cell_num, cont_type, vis_percent)
    %cont_type is who cell_num is contacting, not what cell_num is

    C = get_constants;
    cmap = [1 0 0; .5 .5 .5; 0 0 1; 1 1 0];
    
    c_d = cell_data(cell_num);
    
    if strcmp(cont_type, 'dsgc')
        
        s_conts = detect_vericose_contacts(cell_num, 500, 200, 50000);
        
        dsgc_nums = {[90001 17080 20213 25005 20220], [90002 20125], [20239 20254 20245 20179 20210], [20233 17161 20137 20096]};        
        type_mat = zeros(max(s_conts(1,:)),1);
        for t = 1:length(dsgc_nums)
            type_mat(dsgc_nums{t}) = t;
        end
        
        s_conts(1,:) = type_mat(s_conts(1,:));
        s_conts = s_conts(:,s_conts(1,:) ~= 0);
%         detect_vericose_contacts(sac_nums{l}(s), 500, 200, 50000)
    else
        
        angle2type = zeros(360,1);
        angle2type([1:45, 316:360]) = 1;
        angle2type(46:135) = 2;
        angle2type(136:225) = 3;
        angle2type(226:315) = 4;
        
        unit_circle = [sin((1:360)'/180*pi), cos((1:360)'/180*pi)];
        
        
        dsgc_conts = c_d.contacts;
        cont_nums = unique(dsgc_conts(1,:));
        for n = 1:length(cont_nums)
            if ~any(cont_nums(n) == C.type.(cont_type))
                cont_nums(n) = 0;
            end
        end
        cont_nums(cont_nums==0) = [];
        s_conts = zeros(size(dsgc_conts));
        cont_k = 0;
        
        for n = 1:length(cont_nums)
            sac_d = cell_data(cont_nums(n));
            sac_mid = sac_d.get_midpoint(true);
            
            sub_s_conts = detect_vericose_contacts(cont_nums(n), 500, 200, 50000);
            sub_s_conts = sub_s_conts(:,sub_s_conts(1,:) == cell_num);
            
            cont_xy = [sub_s_conts(4,:)-sac_mid(2); sub_s_conts(5,:)-sac_mid(3)];
            angles = atan2(cont_xy(2,:), cont_xy(1,:));
            angles(angles<0) = angles(angles<0) + 2*pi;
            dir_class = angle2type(ceil(angles/2/pi*360));
            
            s_conts(:, cont_k + (1:length(dir_class))) = [dir_class'; sub_s_conts(2:end,:)];
            cont_k = cont_k + length(dir_class);
            
            
        end
        s_conts = s_conts(:,1:cont_k);
        
    end
    
    figure; plot_cells(cell_num, 1, .01, [0 0 0]);
    
    for n = 1:4
        is_me = s_conts(1,:) == n;
        is_me = is_me & rand(size(is_me)) < vis_percent;
%         sum(is_me)
        scatter(s_conts(4,is_me), s_conts(5,is_me), 50, 'markerFaceColor', cmap(n,:), 'markerEdgeColor', [0 0 0]);
    end
end
        