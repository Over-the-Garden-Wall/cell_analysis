
C = get_constants;


ref_nums = C.type.j(1);
sac_nums = C.type.sure_off_sac;
bins = C.j_bins;
bins = -200:300:100
% bins = 30:60:90

num_bins = length(bins);

% num_angles = 4;
% angle_bins = (0:num_angles)/num_angles*pi;
% angle_bins(end) = angle_bins(end) + .0001;

num_angles = 3;
angle_bins = [0 1/3 2/3 1]*pi;
angle_bins(end) = angle_bins(end) + .0001;


min_thresh = 10000;


% num_angles = length(angle_bins)-1;


load('./j_synapses.mat');
syn_conns(3,:) = 1;
% syn_conns(6,:) = [];


[total_cont, total_denom] = connectivityXdistXangle(ref_nums, sac_nums, bins, angle_bins, syn_conns);
% [total_cont, total_denom] = connectivityXdistXangle(ref_nums, sac_nums, bins, angle_bins, []);
% 




data_mean = zeros(num_bins-1, num_angles);
data_ste = zeros(num_bins-1, num_angles);

for n = 2:num_bins
    for k = 1:num_angles
        is_valid = total_denom(k,n,:) >= min_thresh;
        data_mean(n-1,k) = mean(total_cont(k,n,is_valid)./total_denom(k,n,is_valid))*100;
        data_ste(n-1,k) = std(total_cont(k,n,is_valid)./total_denom(k,n,is_valid))/sqrt(sum(is_valid)-1)*100;
    end
end





figure; hold on;
% errorbar(theta_bins(1:end-1) + bin_size/2, mean_data, ste_data, 'LineWidth', 2);
errorbar(1:num_angles, data_mean, data_ste, 'LineWidth', 2, 'LineStyle', 'none');
bar(1:num_angles, data_mean)
set(gca, 'XTick', 1:3);
set(gca, 'XTickLabel', {'Antagonistic', 'Orthogonal', 'Supporting'});

prep_figure(gcf,gca, 'ylabel', 'Synaptic Contact Density (Synapses / J-voxel)');
ylabel('Synaptic Contact Density (Synapses / J-voxel)','FontSize', 20);

