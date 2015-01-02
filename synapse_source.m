C = get_constants;

cn = C.type.sure_off_sac;

num_sacs = length(cn);
cell_mids = zeros(num_sacs, 3);

j_or_minij = 'minij';

jnum = C.type.(j_or_minij)(1);
proj_plane = 1;
p = 1:3;
p(proj_plane) = [];


for n = 1:num_sacs
    c = cn(n);
    c_d = cell_data(c);
    cell_mids(n,:) = c_d.get_midpoint(true);
end
    
cell_num2ind = sparse(cn,ones(1,num_sacs),1:num_sacs);

if strcmp(j_or_minij,'j')
    load('j_synapses.mat');
    conns = syn_conns;
else
    load('conns_synapses.mat');
end
conns = double(conns);

figure; hold all

plot_cells(jnum, proj_plane,.01,.85*[1 1 1]);

scatter(cell_mids(:,p(1)), cell_mids(:,p(2)), 'o');
% scatter(cell_mids(:,2), cell_mids(:,3), 'o');

is_valid = (conns(1,:) == jnum) | (conns(2,:) == jnum);
conns = conns(:,is_valid);
conns(1:2,conns(2,:) == jnum) = conns([2 1],conns(2,:) == jnum);

synapsing_list = unique(conns(2,:));
syn_inds = cell_num2ind(synapsing_list);
scatter(cell_mids(syn_inds,p(1)), cell_mids(syn_inds,p(2)), '*');
% 
for k = 1:size(conns,2)
    plot([conns(p(1)+3,k); cell_mids(cell_num2ind(conns(2,k)),p(1))], [conns(p(2)+3,k), cell_mids(cell_num2ind(conns(2,k)),p(2))] , 'color', [1 0 1]);
end
scatter(conns(p(1)+3,:), conns(p(2)+3,:))

path_segs = [conns(p(1)+3,:)'-cell_mids(cell_num2ind(conns(2,:)'),p(1)), conns(p(2)+3,k)' - cell_mids(cell_num2ind(conns(2,:)'),p(2))];

figure; hold on
for k = 1:size(path_segs,1)
    plot([0; path_segs(k,1)], [0; path_segs(k,2)]) 
end

% plot(path_segs(:,2), path_segs(:,1));

theta = atan2(path_segs(:,2), path_segs(:,1));
theta = theta+pi/2;
theta(theta > pi) = theta(theta > pi) - 2*pi;
num_bins = 8;
theta_bins = ((1:num_bins)-.5)*2*pi/num_bins - pi;
theta_hist = hist(theta ,theta_bins);

figure; polar([theta_bins theta_bins(1)], [theta_hist theta_hist(1)]);