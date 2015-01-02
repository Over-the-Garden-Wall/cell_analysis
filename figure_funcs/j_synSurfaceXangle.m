C = get_constants
load('/net/omicfs/home/matthew/stratification/j_synapses.mat')

for c = 1:size(syn_conns,2);
    cell_dat = cell_data(syn_conns(2,c));
    cell_mid = cell_dat.get_midpoint(true);
    syn_conns(5,c) = syn_conns(5,c) - cell_mid(2);
    syn_conns(6,c) = syn_conns(6,c) - cell_mid(3);
end

r = sqrt(syn_conns(5,:).^2 + syn_conns(6,:).^2);
x = syn_conns(5,:)./r;
y = syn_conns(6,:)./r;

cell_dat = cell_data(C.type.j(1));
cell_dat.dist_axis

dot_prod = x*cell_dat.dist_axis(1) + y*cell_dat.dist_axis(2);

theta = acos(dot_prod);

contact = syn_conns(3,:);

figure; scatter(contact, theta)

num_bins = 3; 

theta_bins = (0:num_bins)/num_bins * pi;

bin_size = theta_bins(2)-theta_bins(1);

my_bin = ceil(theta/bin_size);

mean_data = zeros(num_bins,1); ste_data = zeros(num_bins,1);

for c = 1:num_bins;
    mean_data(c) = mean(contact(my_bin==c));
    ste_data(c) = std(contact(my_bin==c));
    ste_data(c) = ste_data(c) / sqrt(sum(my_bin==c)-1);
end

figure; hold on;
% errorbar(theta_bins(1:end-1) + bin_size/2, mean_data, ste_data, 'LineWidth', 2);
errorbar(1:num_bins, mean_data, ste_data, 'LineWidth', 2, 'LineStyle', 'none');
bar(1:num_bins, mean_data)
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Antagonistic', 'Orthogonal', 'Supporting'});

prep_figure(gcf,gca, 'ylabel', 'Mean Synaptic Area (edges)');


