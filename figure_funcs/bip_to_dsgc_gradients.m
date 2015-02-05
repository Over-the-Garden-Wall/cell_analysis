C = get_constants;

close all

% oocns = C.type.oodsgc;

oocns = [[90001 25005 20213 20220 17080], [90002], [17161 20233], [20239 20210 20254 20245]];
oocn_dir = [1 1 1 1 1 2 3 3 4 4 4 4];
colors = [1 0 0; 0 .7 0; .8 .8 0; 1 0 1];

num_dsgc = length(oocns);

% oodsgc_axes = zeros(length(oocns),7);
% for n = 1:length(oocns)
%     [mean_axis, pref_axes] = dsgc_preferred_direction(oocns(n));
%     oodsgc_axes(n,:) = [oocns(n), mean_axis', pref_axes(:,1)', pref_axes(:,2)'];    
% end
% save('~/data/stratification/oodsgc_pref_axes.mat', 'oodsgc_axes');

load('~/data/stratification/oodsgc_pref_axes.mat');


types = {{'t2', 't3a', 't3b', 't4'}, {'t5w', 't5h', 't5l'}};

num_bps = cell(2,1);
bp_locs = cell(2,1);
num_types = zeros(2,1);

type_contact = cell(2,1);

plot_x = cell(num_dsgc, sum(num_types));
plot_y = cell(num_dsgc, sum(num_types));

for l = 1:2
    
    num_types(l) = length(types{l});
    num_bps{l} = zeros(num_types(l), 1);
    bp_locs{l} = cell(num_types(l),1);
    
    type_contact{l} = cell(num_types(l),1);
    
    
    for t = 1:num_types(l);
        bp_nums = C.type.(types{l}{t});
        num_bps{l}(t) = length(bp_nums);
        bp_locs{l}{t} = zeros(num_bps{l}(t), 2);
        
        type_contact{l}{t} = zeros(num_bps{l}(t), 1);
        
        for b = 1:num_bps{l}(t)
            c_d = cell_data(bp_nums(b));
            cmid = c_d.get_midpoint(false);
            bp_locs{l}{t}(b,:) = cmid(2:3);
        end
    end
end
            

layer_displacement = zeros(num_dsgc,2);

dtext = {'ventral axis', 'left-right axis'};
daxis = [C.ventral_axis', [0 1; -1 0]*C.ventral_axis'];

       
for dn = 1:num_dsgc
    
    d = oocns(dn);
    
    c_d = cell_data(d);
    
    
    d_p = c_d.get_surface;
    
    p_depth = C.f(d_p(:,1));
    
    p{1} = d_p(p_depth > 10 & p_depth < 50, 2:3);
    p{2} = d_p(p_depth > 50 & p_depth < 80, 2:3);
    
    layer_mid = zeros(2,2);
    
    all_conts = double(c_d.contacts);
    
    
    hulls = cell(2, 1);

    
    for l = 1:2
        

        hull_inds = convhull(p{l}(:,1),p{l}(:,2));
        hulls{l} = [p{l}(hull_inds,1) p{l}(hull_inds,2)];
        
        
        layer_mid(l,:) = mean(p{l});
        
        
        for t = 1:num_types(l)
            
            
            for b = 1:num_bps{l}(t)
                is_me = all_conts(1,:) == C.type.(types{l}{t})(b);
                if any(is_me)
                    type_contact{l}{t}(b) = sum(all_conts(2,is_me));
                else
                    type_contact{l}{t}(b) = 0;
                end
            end
        end
    end
    
    figure('Name', num2str(d), 'Position', [0 0 1200 800]);
    
    x = sum(num_types);        
    y = 1;
    
%     [dummy, pref_axes] = dsgc_preferred_direction(d);
    
    for l = 1:2
        
        [dummy, my_ind] = find(oodsgc_axes(:,1) == d, 1, 'first');        
        daxis = oodsgc_axes(my_ind,l*2 + 2 + (0:1));
        
        for t = 1:num_types(l)
        
            num_vox_in = zeros(num_bps{l}(t),1);
            for b = 1:num_bps{l}(t)
                c_d = cell_data(C.type.(types{l}{t})(b));
                
                h = cell(1,2);
                [h{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
                
                num_vox_in(b) = sum(inpolygon(p{l}(:,1), p{l}(:,2), h{1}, h{2}));
            end
%                 is_in = inpolygon(bp_locs{l}{t}(:,1), bp_locs{l}{t}(:,2), hulls{l}(:,1), hulls{l}(:,2));
            
            
            
            centered_locs = [bp_locs{l}{t}(:,1) - layer_mid(l,1), bp_locs{l}{t}(:,2) - layer_mid(l,2)];
            
            
%             for direction = 1:2
                subplot(y,x, (l-1)*num_types(1) + t); hold on
                
                if sum(num_vox_in) > 0
                    d_along_axis = centered_locs(:,1).*daxis(1) + centered_locs(:,2).*daxis(2);
            
                    x_dat = d_along_axis(num_vox_in>0);
                    y_dat = type_contact{l}{t}(num_vox_in>0) ./ num_vox_in(num_vox_in>0);
                    
                    plot_x{dn, (l-1)*num_types(1) + t} = x_dat;
                    plot_y{dn, (l-1)*num_types(1) + t} = y_dat;
                    
                    [bs, dummy, dummy, dummy, stats] = regress(y_dat,[x_dat ones(size(x_dat))]);

%                     ps{n} = ['p = ' num2str(stats(3))];

                    
%                     bs = regress(y_dat-mean(y_dat),x_dat-mean(x_dat));
%                     
%                     pval = ;
%                     
%                     intercpt = mean(y_dat) - mean(x_dat)*bs;
                    
                    scatter(x_dat, y_dat);
                    plot(x_dat, x_dat*bs(1) + bs(2));
                    
                    title([types{l}{t} ', p=' num2str(stats(3))]);
                    
                else
                    title(['no data, type: ', types{l}{t}]);
                end
%             end
        end
    end
    
    saveas(gcf, ['~/data/stratification/images/bip2dsgc' num2str(d) '.png']);
%     close all


                
            
    
    
end

close all

x_dat = cell(max(oocn_dir),sum(num_types));
y_dat = cell(max(oocn_dir),sum(num_types));
    
for dn = 1:num_dsgc
    for t = 1:sum(num_types)
        x_dat{oocn_dir(dn), t} = [x_dat{oocn_dir(dn), t}; plot_x{dn,t}];
        y_dat{oocn_dir(dn), t} = [y_dat{oocn_dir(dn), t}; plot_y{dn,t}];
    end
end

type_labels = {'BC2', 'BC3a', 'BC3b', 'BC4', 'BC5w', 'BC5o', 'BC5i'};

bar_plot = zeros([size(x_dat,2), 2]);
bar_error = zeros([size(x_dat,2), 2]);


for t = 1:size(x_dat,2)
    figure; hold all; title(type_labels{t});

    h = zeros(size(x_dat,1)+1,1);
    ps = cell(size(x_dat,1)+1,1);
    tot_x = [];
    tot_y = [];
    for n = 1:size(x_dat,1)
        
        
        scatter(x_dat{n,t}, y_dat{n,t}, '*', 'markerEdgeColor', colors(n,:));        
        tot_x = [tot_x; x_dat{n,t}];
        tot_y = [tot_y; y_dat{n,t}];
        
        [bs, dummy, dummy, dummy, stats] = regress(y_dat{n,t},[x_dat{n,t} ones(size(x_dat{n,t}))]);

        ps{n} = ['p = ' num2str(stats(3))];

        intercpt = mean(y_dat{n,t}) - mean(x_dat{n,t})*bs;

        h(n) = plot(x_dat{n,t}, x_dat{n,t}*bs(1) + bs(2), 'lineWidth', 2, 'Color', colors(n,:));
    end
    
    
        bar_plot(t, 1) = mean(tot_y(tot_x<0));
        bar_plot(t, 2) = mean(tot_y(tot_x>0));
        
        bar_error(t, 1) = std(tot_y(tot_x<0)) / sqrt(sum(tot_x<0));
        bar_error(t, 2) = std(tot_y(tot_x>0)) / sqrt(sum(tot_x>0));
        
    
    [bs, dummy, dummy, dummy, stats] = regress(tot_y,[tot_x ones(size(tot_x))]);

    ps{size(x_dat,1)+1} = ['p = ' num2str(stats(3))];
        

    h(size(x_dat,1)+1) = plot(tot_x, tot_x*bs(1) + bs(2), 'lineWidth', 2, 'Color', [0 0 0]);
    
    legend(h, ps);
    
    saveas(gcf, ['~/data/stratification/images/bip2dsgc_' type_labels{t} '.png']);
    
end


    figure;
    errorbar(bar_plot, bar_error);
% end
