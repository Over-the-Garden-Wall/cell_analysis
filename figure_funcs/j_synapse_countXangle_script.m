load('/net/omicfs/home/matthew/stratification/j_synapses.mat');

C = get_constants;
conns = syn_conns;
num_conns = size(conns,2);

% load(C.trans_loc);
% for n = 1:num_conns
%     conns(6:8,n) = apply_transform(T,conns(6:8,n)');
% end


js = unique(conns(1,:));
sacs = C.type.sure_off_sac;


theta_bins = (0:3)/3*pi;
num_bins = length(theta_bins)-1;


j_axes = zeros(length(js),2);
for n = 1:length(js)
    cell_dat = cell_data(js(n));
    j_axes(n,:) = cell_dat.dist_axis;
end

jid2ind = sparse(js, ones(size(js)), 1:length(js));


sac_mids = zeros(length(sacs),3);
for n = 1:length(sacs)
    cell_dat = cell_data(sacs(n));
    sac_mids(n,:) = cell_dat.get_midpoint(true);
end
sac_mids(:,1) = [];

sacid2ind = sparse(sacs, ones(size(sacs)), 1:length(sacs));

conn_rel_loc = zeros(num_conns,2);

conn_dist = zeros(num_conns,1);
conn_theta = zeros(num_conns,1);


for n = 1:num_conns
    sac_id = conns(2,n);
    jind = jid2ind(conns(1,n));
    
    for d = 1:2;
        conn_rel_loc(n,d) = conns(4+d,n) - sac_mids(sacid2ind(sac_id),d);
    end
    conn_dist(n) = sqrt(sum(conn_rel_loc(n,:).^2));
    
    loc_vec = conn_rel_loc(n,:)/conn_dist(n);
    loc_dot = sum(loc_vec.*j_axes(jind,:));
    
    conn_theta(n) = abs(acos(loc_dot));
end

syn_bin_hist = zeros(num_bins,1);
for n = 1:num_bins
    is_this_angle = (conn_theta>=theta_bins(n)) & (conn_theta < theta_bins(n+1));
    syn_bin_hist(n) = sum(is_this_angle);
end


min_dist = min(conn_dist);
max_dist = max(conn_dist);




denom_bin_hist = zeros(num_bins,1);
    
for n = 1:num_conns
    rel_loc_to_sac = [conns(5,n)-sac_mids(:,1), conns(6,n)-sac_mids(:,2)];
    rel_dist = sqrt(sum(rel_loc_to_sac.^2,2));
    dist_is_good = rel_dist >= min_dist & rel_dist <= max_dist;
    rel_vec = [rel_loc_to_sac(:,1)./rel_dist, rel_loc_to_sac(:,2)./rel_dist];
    my_j = jid2ind(conns(1,n));
    rel_dot = rel_vec(:,1)*j_axes(my_j,1) + rel_vec(:,2)*j_axes(my_j,2);
    rel_theta = abs(acos(rel_dot));
    rel_theta = rel_theta(dist_is_good);
    
    
    for k = 1:num_bins
        is_this_angle = (rel_theta>=theta_bins(k)) & (rel_theta < theta_bins(k+1));
        denom_bin_hist(k) = denom_bin_hist(k)+sum(is_this_angle);
    end
end
    
bin_ratio = syn_bin_hist./denom_bin_hist;

bar(bin_ratio/sum(bin_ratio));
set(gca, 'XTick', 1:3);
set(gca, 'YTick', .1:.1:1);
set(gca, 'YLim', [0 .5]);

set(gca, 'XTickLabel', {'Antagonistic', 'Orthogonal', 'Supporting'});

prep_figure(gcf,gca);
ylabel('Adjusted Probability','FontSize', 20);





