C = get_constants;

cn = C.type.sure_off_sac;
% types = {'sure_off_sac', 't1', 't2', 't3a', 't3b', 't4'};
types = {'sure_off_sac', 't1', 't2'};
typelabels = {'Off SAC', 'BC1', 'BC2'};
% types = {'sure_off_sac', 't2'};
num_types = length(types);

num_sacs = length(cn);
cell_mids = zeros(num_sacs, 3);

j_or_minij = 'j';
%minijs: 1 2 5
%js: 1 2
jnum = C.type.(j_or_minij)(1);


for n = 1:num_sacs
    c = cn(n);
    c_d = cell_data(c);
    cell_mids(n,:) = c_d.get_midpoint(true);
end
c_d = cell_data(jnum);

p = c_d.get_surface;
ax = c_d.dist_axis;
j_mid = c_d.get_midpoint(true);


theta = atan2(ax(1),ax(2));
Q = [1 0 0; 0 cos(theta) -sin(theta); 0 sin(theta) cos(theta)];
p = p*Q';
rot_jmid = j_mid*Q';

cell_num2ind = sparse(cn,ones(1,num_sacs),1:num_sacs);

if strcmp(j_or_minij,'j')
    load('j_synapses.mat');
    sac_conns = syn_conns;
else
    load('conns_synapses.mat');
    sac_conns = conns;
end

sac_conns = double(sac_conns);

is_valid = (sac_conns(1,:) == jnum) | (sac_conns(2,:) == jnum);
sac_conns = sac_conns(:,is_valid);
sac_conns(1:2,sac_conns(2,:) == jnum) = sac_conns([2 1],sac_conns(2,:) == jnum);
sac_conns(4:6,:) = Q*sac_conns(4:6,:);
bip_conns = double(c_d.contacts);
bip_conns(3:5,:) = Q*bip_conns(3:5,:);

conns = cell(length(types),1);
conns{1} = sac_conns(2:end,:);

for k = 2:num_types
    type_num = false(1,size(bip_conns,2));
    for t = 1:size(bip_conns,2)
        if any(bip_conns(1,t)==C.type.(types{k}))
            type_num(t) = true;
        end
    end
    conns{k} = bip_conns(:,type_num);
end


num_p = size(p,1);


% plot_cells(jnum, proj_plane,.01,.85*[1 1 1]);

% scatter(cell_mids(:,p(1)), cell_mids(:,p(2)), 'o');
% scatter(cell_mids(:,2), cell_mids(:,3), 'o');

% synapsing_list = unique(sac_conns(2,:));
% syn_inds = cell_num2ind(synapsing_list);
% scatter(cell_mids(syn_inds,p(1)), cell_mids(syn_inds,p(2)), '*');
% 
% for k = 1:size(conns,2)
%     plot([conns(6,k); cell_mids(cell_num2ind(conns(2,k)),2)], [conns(5,k), cell_mids(cell_num2ind(conns(2,k)),3)] , 'y');
% end

figure; hold all
scatter(p(1:100:num_p,2), p(1:100:num_p,3), 1, 'markerEdgeColor', [1 1 1]*.85);

h = zeros(num_types,1);
for k = 1:num_types    
    h(k) = scatter(conns{k}(4,:), conns{k}(5,:), '*');        
end
legend(h, typelabels);


figure; hold all
scatter(p(1:100:num_p,3), p(1:100:num_p,1), 1, 'markerEdgeColor', [1 1 1]*.85);

h = zeros(num_types,1);
for k = 1:num_types    
    h(k) = scatter(conns{k}(5,:), conns{k}(3,:), '*');        
end
legend(h, typelabels);

contact_sum = zeros(num_types-1,1);
for k = 2:num_types
    contact_sum(k-1) = sum(conns{k}(2,:))/length(unique(conns{k}(1,:)));
end
figure; bar(contact_sum)

bin_size = 10000;
bins = -bin_size:bin_size:max(bip_conns(5,:)-rot_jmid(3));
dist_hists = zeros(length(bins), num_types);
dist_cont_tot = zeros(length(bins), num_types);
for k = 1:num_types
    dists = conns{k}(5,:) - rot_jmid(3);
    
    dist_hists(:,k) = hist(dists,bins);
    
    for n = 1:length(bins)
        is_valid = dists> bins(n)-bin_size/2 & dists< bins(n)+bin_size/2;
        dist_cont_tot(n,k) = sum(conns{k}(2,is_valid));
    end
end

% figure; hold all
% for k = 1:num_types
%     plot(bins, dist_hists(:,k), 'LineWidth', 2);
% end
%     legend(typelabels);
    
    
    

figure; hold all
for k = 1:num_types
    plot(bins, dist_cont_tot(:,k), 'LineWidth', 2);
end
    legend(typelabels);
