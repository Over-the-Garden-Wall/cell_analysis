function plot_sac_contacts(cell_num, size_thresh)

    C = get_constants;
    
%     sk = 8;

    cell_dat = cell_data(cell_num);
    
    conts = double(cell_dat.contacts);
    
    warning('hack')
    load('./j_synapses');
    conns = syn_conns;
    conns(:,conns(1,:)~=cell_num) = [];
    conts = conns(2:end,:);
    
    
    
    is_sac = false(size(conts,2),1);
    for k = 1:size(conts,2)
        is_sac(k) = any(conts(1,k)==C.type.sure_off_sac);
    end
    
    conts = conts(:,is_sac);
    
    num_conts = size(conts,2);
    
    arrow_size = conts(2,:)';
    arrow_size(arrow_size < size_thresh) = 0;
    
    arrow_loc = conts(4:5,:)';
    
    arrow_dir = zeros(num_conts,2);
    sac_nums = conts(1,:);
    
    for k = 1:num_conts
        sac_data = cell_data(sac_nums(k));
        sac_mid = sac_data.get_midpoint(true);
        arrow_dir(k,:) = sac_mid(2:3) - arrow_loc(k,:);
        arrow_dir(k,:) = arrow_dir(k,:) / sqrt(sum(arrow_dir(k,:).^2));
    end
    
    
    figure;
    plot_cells(cell_num, 1, .01, [0 0 0]);
    hold on
%     plot_cells(C.type.sure_off_sac(sk), 1, .01, [1 0 0]);
    
    
    
    quiver(arrow_loc(:,1), arrow_loc(:,2), arrow_size.*arrow_dir(:,1), arrow_size.*arrow_dir(:,2)); 
    
    title(num2str(cell_num));
    
    figure; hist(atan2(arrow_dir(:,2),arrow_dir(:,1))); 
    
    
end
    
    
   
    