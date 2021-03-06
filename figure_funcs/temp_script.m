C = get_constants;

c{1} = C.type.sure_off_sac;
c{2} = C.type.on_sac;

for n = 1:length(c);
    cinfo = contact_info(c{n},c{n}, false);
    
    pairs = unique(cinfo.cell_ids, 'rows');
    soma_locs = [cinfo.soma_loc(pairs(:,1),:), cinfo.soma_loc(pairs(:,2),:)];
    p_dist = sqrt(sum((soma_locs(:,1:3)-soma_locs(:,4:6)).^2 , 2));
    
    p_total = zeros(size(pairs,1),1);
    for p = 1:size(pairs,1);
        is_me = cinfo.cell_ids(:,1) == pairs(p,1) & cinfo.cell_ids(:,2) == pairs(p,2);
        p_total(p) = sum(cinfo.contact_size(is_me));
    end
    
    [p_whist x{n}] = weighted_hist(p_dist, p_total, 15);
    p_hist = hist(p_dist, x{n});
    
    p_data{n} = p_whist ./ p_hist;
end

figure; hold all
for n = 1:length(p_data);
    plot(x{n}, p_data{n}, 'lineWidth', 2);
end