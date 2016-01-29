%bc-gc contact script.


type_is_predefined = false;
output_is_coverage = true;
is_on_dsgc = false;
use_face_only = false;


% dsgcs = C.type.on_sac;

bc_to_dsgc_contxdistxoverlap

cns = [];
for n = 1:length(BC_types)
    cns{n} = C.type.(BC_types{n});
end

a = peters_rule_prediction(dsgcs, cns, [0 70]);

densities = type_contact./type_overlap;

mean_dens = zeros(1,size(densities,2));
std_dens = zeros(1,size(densities,2));

densities(isinf(densities)) = NaN;

for n = 1:size(densities,2)
    mean_dens(n) = mean(densities(~isnan(densities(:,n)),n));
    std_dens(n) = std(densities(~isnan(densities(:,n)),n));    
end
ste_dens = std_dens / sqrt(size(densities,1)-1);
figure; 

h = error_dot_plot(100*[mean_dens' a], 100*[ste_dens' zeros(10,1)], BC_type_lbls);

legend(h, {'Observed Contact', 'Expected Contact'});
prep_figure(gcf,gca)
ylabel('Percent DSGC Coverage', 'fontSize', 12/.45);

pr_pred = zeros(size(densities));
for n = 1:size(densities,1)
    pr_pred(n,:) = peters_rule_prediction(dsgcs(n), cns, [0 70])';
end

figure; 

h = error_dot_plot(100*(densities - pr_pred)', zeros(size(densities')), BC_type_lbls);

lgnd = cell(1, size(densities,1));
for k = 1:length(lgnd)
    lgnd{k} = num2str(dsgcs(k));
end

legend(lgnd);
prep_figure(gcf,gca);
ylabel('Deviation from PR prediction (%coverage)', 'fontSize', 12/.45);

