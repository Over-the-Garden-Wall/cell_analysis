%test excitatory DS hypothesis

%the hypothesis is that bipolar inputs to the ooDSGCs are DS as a result of
%directionally specific SAC input to BC axon terminals.
close all
C = get_constants;

SAC_output_zone = 100000;
terminal_distance_threshold = 1000;
num_angle_bins = 12;

exclude_touching_sacs = true;

angle_bins = (1:num_angle_bins)/num_angle_bins * 2 * pi;

CELL_MAX = 100000;

BC_types = {'t2', 't3a', 't3b', 't4', 't5w', 't5h', 't5l'};
% BC_types = {'t3a'};

% directions = [atan2(0,1); atan2(1,0); atan2(0, -1); atan2(-1, 0)];
directions = [atan2(1,0); atan2(0,1); atan2(-1, 0); atan2(0, -1)];

ooDSGCs = [17161, 20233, ... %down PD
    20239, 20254, 20210, 20245, 20179, ... right PD
    90002, ... %up PD
    25005, 90001, 20213, 17080, 20220]; %left PD

ooDSGC_to_direction = [3, 3, 2, 2, 2, 2, 2, 1, 4, 4, 4, 4, 4];

is_ooDSGC = zeros(CELL_MAX,1);
is_ooDSGC(ooDSGCs) = 1:length(ooDSGCs);

SACs = [C.type.sure_off_sac, C.type.on_sac];
sac_ind = zeros(CELL_MAX,1);
sac_ind(SACs) = 1:length(SACs);


sac_mids = zeros(length(SACs),2);

%get sac midpoints and map sacs
for n = 1:length(SACs);
    c_d = cell_data(SACs(n));
    t = c_d.get_midpoint;
    sac_mids(n,:) = t(2:3);
end


bin_count = zeros(length(angle_bins), length(BC_types), length(ooDSGCs));
bin_contact = zeros(length(angle_bins), length(BC_types), length(ooDSGCs));

dsgc_bc_count = zeros(length(ooDSGCs),length(BC_types));
dsgc_bc_last = zeros(length(ooDSGCs),length(BC_types));
        
for type_k = 1:length(BC_types)
    
    cns = C.type.(BC_types{type_k});
    
    for ck = 1:length(cns)
        c = cns(ck);
        c_d = cell_data(c);
        
        my_contacts = double(c_d.contacts);
        
        my_oodsgc_contacts = my_contacts(:,is_ooDSGC(my_contacts(1,:))>0);
        my_sac_contacts = my_contacts(:,sac_ind(my_contacts(1,:)) > 0);
        
        valid_sac_cont = false(1,size(my_sac_contacts,2));
        
        for sck = 1:size(my_sac_contacts,2)
            midp = sac_mids(sac_ind(my_sac_contacts(1, sck)),:);
            d = sqrt( sum( (midp - my_sac_contacts(4:5, sck)').^2 ) );
            valid_sac_cont(sck) = d > SAC_output_zone;
        end
            
        my_sac_contacts = my_sac_contacts(:, valid_sac_cont);
        
        oodsgc_linked_contacts = zeros(size(my_oodsgc_contacts,2), 1);
        
        for n = 1:size(my_oodsgc_contacts,2)
            ds = sqrt(sum((my_sac_contacts(3:5, :) - ...
                (my_oodsgc_contacts(3:5, n) * ...
                ones(1,size(my_sac_contacts,2)))).^2));
            [min_d, min_ind] = min(ds);
            if min_d <= terminal_distance_threshold
                oodsgc_linked_contacts(n) = min_ind;
                
                if exclude_touching_sacs
                   c_d = cell_data( my_oodsgc_contacts(1,n) );
                   from_oodsgc_conts = double(c_d.contacts);
                   sac_num = my_sac_contacts(1, oodsgc_linked_contacts(n));
                   poss_excluding_contact = from_oodsgc_conts(1,:) == sac_num;
                   ds = sqrt(sum((from_oodsgc_conts(3:5, poss_excluding_contact) - ...
                        (my_oodsgc_contacts(3:5, n) * ...
                        ones(1,sum(poss_excluding_contact)))).^2));

                    if any(ds < terminal_distance_threshold * 2)
                       oodsgc_linked_contacts(n) = 0;
                   end
                    
                end
                
            end
        end
        
        for n = 1:length(oodsgc_linked_contacts)
            if oodsgc_linked_contacts(n) > 0
                sac_num = my_sac_contacts(1, oodsgc_linked_contacts(n));
                DSGC_num = is_ooDSGC(my_oodsgc_contacts(1,n));
                
                oodsgc_dir = directions(ooDSGC_to_direction(DSGC_num));
                
                
                sac_m = sac_mids(sac_ind(sac_num),:);
                cont_loc = my_sac_contacts(4:5, oodsgc_linked_contacts(n))';
                
                sac_angle = atan2(cont_loc(2) - sac_m(2), cont_loc(1) - sac_m(1));
                
                angle_diff = sac_angle - oodsgc_dir;
                if angle_diff < 0
                    angle_diff = angle_diff + 2*pi;
                end
                
                bin_num = ceil(angle_diff / 2 / pi * num_angle_bins);
                
                bin_count(bin_num, type_k, DSGC_num) = bin_count(bin_num, type_k, DSGC_num) + 1;
                bin_contact(bin_num, type_k, DSGC_num) = bin_contact(bin_num, type_k, DSGC_num) + ...
                    my_oodsgc_contacts(2,n);
                
                if dsgc_bc_last(DSGC_num, type_k) ~= c                    
                    dsgc_bc_count(DSGC_num, type_k) = dsgc_bc_count(DSGC_num, type_k) + 1;
                    dsgc_bc_last(DSGC_num, type_k) = c;
                end
                
                if type_k == 2 && (bin_num <= 3 || bin_num >= 11)
                    fprintf('%d ', [sac_num, ooDSGCs(DSGC_num), c, ...
                        my_oodsgc_contacts(2:5,n)']);
                    fprintf('\n');
                end
                        
                
            end
        end
        
        
    end
    
    
end


%
ooDSGCs_bin_count = zeros(length(angle_bins), 1);
ooDSGCs_bin_contact = zeros(length(angle_bins), 1);

for k = 1:length(ooDSGCs)
   c_d = cell_data(ooDSGCs(k));
   my_dir = directions(ooDSGC_to_direction(k));
   
   my_conts = double(c_d.contacts);
   my_conts = my_conts(:,sac_ind(my_conts(1,:))>0);
   
   for n = 1:size(my_conts, 2)
       sac_m = sac_mids(sac_ind(my_conts(1,n)),:);
       cont_loc = my_conts(4:5, n)';

       sac_angle = atan2(cont_loc(2) - sac_m(2), cont_loc(1) - sac_m(1));

       angle_diff = sac_angle - my_dir;
       if angle_diff < 0
           angle_diff = angle_diff + 2*pi;
       end

       bin_num = ceil(angle_diff / 2 / pi * num_angle_bins);

       ooDSGCs_bin_count(bin_num) = ooDSGCs_bin_count(bin_num) + 1;
       ooDSGCs_bin_contact(bin_num) = ooDSGCs_bin_contact(bin_num) + ...
           my_conts(2,n);
    end
end
%

angles = (0:num_angle_bins)' / num_angle_bins * 2 * pi;

figure; 

polar(angles * ones(1,size(bin_count,2)), ...
    sum(bin_count([1:end 1],:,:),3));

legend(BC_types);

figure;
polar(angles, ...
    ooDSGCs_bin_count([1:end 1]));

angles(end) = [];
angles = angles / pi * 180;
bar_plot = zeros(length(BC_types),3);
for n = 1:length(BC_types)
    bar_plot(n,1) = mean(sum(bin_contact(angles <= 60 | angles >= 300, n,:),3));
    bar_plot(n,2) = ...
        mean(sum(bin_contact((angles < 120 & angles > 60) | ...
        (angles < 300 & angles > 240), n,:),3));
    bar_plot(n,3) = mean(sum(bin_contact(angles <= 240 & angles >= 120, n,:),3));
end

figure; bar(bar_plot);
set(gca, 'XTickLabel', BC_types);
legend({'ND', 'Orthogonal', 'PD'});
    
ylabel('contact area');


error_bar_plot = zeros(length(BC_types), 2);
is_ND = angles <= 60 | angles >= 300;
is_PD = angles <= 240 & angles >= 120;

for n = 1:length(BC_types)
    for k = 1:size(dsgc_bc_count,1)
        bin_contact(:, n, k) = bin_contact(:, n, k)/ dsgc_bc_count(k, n);        
    end
    bin_contact(isnan(bin_contact)) = 0;
    PD_contact = squeeze(sum(bin_contact(is_PD, n, :),1));
    ND_contact = squeeze(sum(bin_contact(is_ND, n, :),1));
    error_bar_plot(n,1) = mean(PD_contact - ND_contact);
    error_bar_plot(n,2) = std(PD_contact - ND_contact) / sqrt(size(bin_contact,3)-1);
end
figure;
errorbar(error_bar_plot(:,1), error_bar_plot(:,2), 'lineWidth', 2); 
hold on; plot([.5, size(error_bar_plot,1)+.5], [0, 0], '--k');
set(gca, 'XTick', 1:size(error_bar_plot,1));
set(gca, 'XTickLabel', BC_types);
ylabel('mean contact difference')
    