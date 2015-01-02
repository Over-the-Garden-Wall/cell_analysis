C = get_constants;


ref_nums = C.type.j;
sac_nums = C.type.sure_off_sac;
bins = C.j_bins;
num_bins = length(bins);

num_angles = 4;
angle_bins = (0:num_angles)/num_angles*pi;
angle_bins(end) = angle_bins(end) + .0001;

min_thresh = 10000;


% num_angles = length(angle_bins)-1;


% load('./j_synapses.mat');
% syn_conns(3,:) = 1;
% % syn_conns(6,:) = [];
% 
% 
% [total_cont, total_denom] = connectivityXdistXangle(ref_nums, sac_nums, bins, angle_bins, syn_conns);
[total_cont, total_denom] = connectivityXdistXangle(ref_nums, sac_nums, bins, angle_bins, []);





data_mean = zeros(num_bins, num_angles);
data_ste = zeros(num_bins, num_angles);

for n = 1:num_bins
    for k = 1:num_angles
        is_valid = total_denom(k,n,:) >= min_thresh;
        data_mean(n,k) = mean(total_cont(k,n,is_valid)./total_denom(k,n,is_valid))*100;
        data_ste(n,k) = std(total_cont(k,n,is_valid)./total_denom(k,n,is_valid))/sqrt(sum(is_valid)-1)*100;
    end
end

        




figure; hold all

leg_entries = cell(num_angles,1);
for k = 1:num_angles;
    errorbar(bins, data_mean(:,k), data_ste(:,k), 'LineWidth', 2);
    leg_entries{k} = [num2str(round(angle_bins(k)/pi*180)) ' to ' num2str(round(angle_bins(k+1)/pi*180))];
end


prep_figure(gcf,gca, 'legend', leg_entries, ...
    'xlabel', 'Distance from soma (microns)', ...
    'ylabel', 'contact (%)');




