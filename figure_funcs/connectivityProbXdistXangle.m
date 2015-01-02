function [total_cont, total_denom, total_denom_sac] = connectivityProbXdistXangle(ref_nums, conn_nums, dist_bins, angle_bins, conns)

C = get_constants;





% ref_nums = C.type.sure_off_sac;
% ref_nums = C.type.minij;
num_refs = length(ref_nums);

% conn_nums = C.type.sure_off_sac;
num_conn_cells = length(conn_nums);



% dist_bins = C.minij!b;
bin_size = dist_bins(2)-dist_bins(1);
num_bins = length(dist_bins);

% angle_bins = (0:num_angles)/num_angles*pi;
% angle_bins(end) = angle_bins(end) + .0001;
num_angles = length(angle_bins)-1;

total_cont = zeros( num_angles, num_bins, num_conn_cells, num_refs);
total_denom = zeros( num_angles, num_bins, num_conn_cells, num_refs);
total_denom_sac = zeros( num_angles, num_bins, num_conn_cells, num_refs);

% load('./j_synapses.mat');

for r = 1:num_refs
    
    
%     figure; cmap = colormap('lines'); plot_cells(ref_nums(r),1,.01,[0 0 0]); hold all;
    
    
    
    cell_dat = cell_data(ref_nums(r));
    p = cell_dat.get_surface;
    
    is_hull = convhull(p(:,2), p(:,3));
    j_hull = p(is_hull, [2 3]);
    
    cell_mid = cell_dat.get_midpoint(true);
    
    tp = [p(:,2)-cell_mid(2), p(:,3)-cell_mid(3)];
        
    d = tp(:,1)*cell_dat.dist_axis(1) + ...
        tp(:,2)*cell_dat.dist_axis(2);
    d = d/1000;        
        
    d = d - dist_bins(1) + bin_size/2;
    d = ceil(d/bin_size);
        
    is_valid = d >= 1 & d <= num_bins;
    d = d(is_valid);
    p = p(is_valid,:);
    
    phi = atan2(cell_dat.dist_axis(2), cell_dat.dist_axis(1));
    
    
    if ~exist('conns', 'var') || isempty(conns)
        conts = double(cell_dat.contacts);
    else
        conts = conns(2:end, conns(1,:) == ref_nums(r));
    end
    
    
    
    conts(6,:) = ((conts(4,:)-cell_mid(2))*cell_dat.dist_axis(1) + ...
            (conts(5,:)-cell_mid(3))*cell_dat.dist_axis(2))/1000;
    conts(6,:) = ceil((conts(6,:) - dist_bins(1) + bin_size/2)/bin_size);
    is_valid = conts(6,:) >= 1 & conts(6,:) <= num_bins;
    conts = conts(:,is_valid);
        
    for s = 1:num_conn_cells
        sac_dat = cell_data(conn_nums(s));
        
        sac_p = sac_dat.get_surface;
        
        sac_p = sac_p(inpolygon(sac_p(:,2), sac_p(:,3), j_hull(:,1), j_hull(:,2)),:);
        
        
        
        sac_mid = sac_dat.get_midpoint(true);
        
        sac_tp = [sac_p(:,2)-cell_mid(2), sac_p(:,3)-cell_mid(3)]; 
        sac_d = sqrt(sac_tp(:,1).^2 + sac_tp(:,2).^2);
        
        sac_tp = [sac_p(:,2)-sac_mid(2), sac_p(:,3)-sac_mid(3)]; 
        sac_td = sqrt(sac_tp(:,1).^2 + sac_tp(:,2).^2);
        sac_theta = acos(sac_tp(:,1)./sac_td*cell_dat.dist_axis(1) + sac_tp(:,2)./sac_td*cell_dat.dist_axis(2));
        
        is_valid = sac_td < 130*1000 & sac_td > 80*1000;
        sac_theta = sac_theta(is_valid);
        sac_d = sac_d(is_valid)/1000;
        
        sac_d = sac_d - dist_bins(1) + bin_size/2;
        sac_d = ceil(sac_d/bin_size);
        
        
        
        tp = [p(:,2)-sac_mid(2), p(:,3)-sac_mid(3)];       
        td = sqrt(tp(:,1).^2 + tp(:,2).^2);
        is_valid = td < 130*1000 & td > 80*1000;

        tp = [tp(is_valid,1)./td(is_valid) tp(is_valid,2)./td(is_valid)];
        
        theta = acos(tp(:,1)*cell_dat.dist_axis(1) + tp(:,2)*cell_dat.dist_axis(2));
        
        
%         theta = atan2(tp(is_valid,2), tp(is_valid,1));
        
%         theta = mod(abs(theta - phi),pi);

        for k = 1:num_angles
            angle_flag = theta >= angle_bins(k) & theta < angle_bins(k+1);
            sac_angle_flag = sac_theta >= angle_bins(k) & sac_theta < angle_bins(k+1);
            
            for n = 1:num_bins
                total_denom(k, n, s, r) = sum(angle_flag & d(is_valid)==n);
                total_denom_sac(k, n, s, r) = sum(sac_angle_flag & sac_d==n);
                
            end
        end
    
        my_conts = conts(:,conts(1,:)==conn_nums(s));
        
        
        norm_conts = [my_conts(4,:) - sac_mid(2); my_conts(5,:) - sac_mid(3)];
        cont_r = sqrt(norm_conts(1,:).^2 + norm_conts(2,:).^2);
        norm_conts(1,:) = norm_conts(1,:)./cont_r;
        norm_conts(2,:) = norm_conts(2,:)./cont_r;
        
        if ~isempty(my_conts)
            cont_theta = acos(norm_conts(1,:)*cell_dat.dist_axis(1) + ...
                (norm_conts(2,:)*cell_dat.dist_axis(2)));
            
            cont_bins = zeros(size(cont_theta));
            for k = 1:num_angles
                angle_flag = cont_theta >= angle_bins(k) & cont_theta < angle_bins(k+1);
                cont_bins(angle_flag) = k;
            end
                
            
            
            for c = 1:size(my_conts,2)
                total_cont(cont_bins(c), my_conts(6,c), s, r) = ...
                    total_cont(cont_bins(c), my_conts(6,c), s, r) + my_conts(2,c);
                
                
%                 plot([sac_mid(2); my_conts(4,c)], [sac_mid(3); my_conts(5,c)], 'LineWidth', 2, 'Color', cmap(cont_bins(c),:));
                
            end
            
            
        end
        
    end
    
    
        
    
end



end




