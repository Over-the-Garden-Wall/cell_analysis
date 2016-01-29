C = get_constants;

close all

% oocns = C.type.oodsgc;

oocns = [[90001 25005 20213 20220 17080], [90002], [17161 20233], [20239 20210 20254 20245], 20179];
oocn_dir = [1 1 1 1 1 2 3 3 4 4 4 4 5];
colors = [1 0 0; 0 .7 0; .8 .8 0; 1 0 1; 0 0 0];

num_dsgc = length(oocns);

% oodsgc_axes = zeros(length(oocns),7);
% for n = 1:length(oocns)
%     [mean_axis, pref_axes] = dsgc_preferred_direction(oocns(n));
%     oodsgc_axes(n,:) = [oocns(n), mean_axis', pref_axes(:,1)', pref_axes(:,2)'];    
% end
% save('~/data/stratification/oodsgc_pref_axes.mat', 'oodsgc_axes');

load('~/data/stratification/oodsgc_pref_axes.mat');


types = {'sure_off_sac', 'on_sac'};

num_types = length(types);

type_contact = cell(num_types,num_dsgc);
type_loc = cell(num_types,num_dsgc);


plot_x = cell(num_dsgc, sum(num_types));
plot_y = cell(num_dsgc, sum(num_types));    


       
for dn = 1:num_dsgc
    
    d = oocns(dn);
    
    c_d = cell_data(d);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
    p{1} = d_p(p_depth > 10 & p_depth < 50, 2:3);
    p{2} = d_p(p_depth > 50 & p_depth < 80, 2:3);
    
%     layer_mid = zeros(2,2);
    
    all_conts = double(c_d.contacts);
    
%     for l = 1:2;
    
        
        for t = 1:num_types
            [dummy, my_ind] = find(oodsgc_axes(:,1) == d, 1, 'first');   
%             daxis = oodsgc_axes(my_ind,l*2 + 2 + (0:1));
            daxis = oodsgc_axes(my_ind,4:5);
            
            layer_mid = mean(p{t});
            
            cns = C.type.(types{t});
            is_my_type = false(1,size(all_conts,2));
            for b = 1:length(cns)
                is_my_type = is_my_type | (all_conts(1,:) == cns(b));                
            end
            my_conts = all_conts(:,is_my_type);
            my_conts(2,:) = min(my_conts(2,:), 1000);
            
            dist = (my_conts(4,:) - layer_mid(1))*daxis(1) + (my_conts(5,:) - layer_mid(2))*daxis(2);
            type_loc{t, dn} = dist';
            type_contact{t,dn} = my_conts(2,:)';
        end
    
%     figure('Name', num2str(d), 'Position', [0 0 1200 800]);
     
            
    
    
end

close all

x_dat = cell(max(oocn_dir),sum(num_types));
y_dat = cell(max(oocn_dir),sum(num_types));
    
for dn = 1:num_dsgc
    for t = 1:sum(num_types)
        x_dat{oocn_dir(dn), t} = [x_dat{oocn_dir(dn), t}; type_loc{t, dn}];
        y_dat{oocn_dir(dn), t} = [y_dat{oocn_dir(dn), t}; type_contact{t, dn}];
    end
end

type_labels = {'OFF SAC', 'ON SAC'};

for t = 1:size(x_dat,2)
    figure; hold all; title(type_labels{t});

    h = zeros(size(x_dat,1)+1,1);
    ps = cell(size(x_dat,1)+1,1);
    tot_x = [];
    tot_y = [];
    for n = 1:size(x_dat,1)
%         scatter(x_dat{n,t}, y_dat{n,t}, '*', 'markerEdgeColor', colors(n,:));        
        tot_x = [tot_x; x_dat{n,t}];
        tot_y = [tot_y; y_dat{n,t}];
        
        [bs, dummy, dummy, dummy, stats] = regress(y_dat{n,t},[x_dat{n,t} ones(size(x_dat{n,t}))]);

        ps{n} = ['p = ' num2str(stats(3))];

        intercpt = mean(y_dat{n,t}) - mean(x_dat{n,t})*bs;

        h(n) = plot(x_dat{n,t}, x_dat{n,t}*bs(1) + bs(2), 'lineWidth', 2, 'Color', colors(n,:));
    end
    
    [bs, dummy, dummy, dummy, stats] = regress(tot_y,[tot_x ones(size(tot_x))]);

    ps{size(x_dat,1)+1} = ['p = ' num2str(stats(3))];
        

    h(size(x_dat,1)+1) = plot(tot_x, tot_x*bs(1) + bs(2), 'lineWidth', 2, 'Color', [0 0 0]);
    
    legend(h, ps);
end
