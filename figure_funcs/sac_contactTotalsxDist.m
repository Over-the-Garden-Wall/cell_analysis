C = get_constants;

types = {'t1', 't2', 't3a', 't3b', 't4'};

sac_nums = C.type.sure_off_sac;

num_sacs = length(sac_nums);
num_types = length(types);

bins = C.sac_bins;
bin_size = bins(2)-bins(1);
num_bins = length(bins);

% figure; main_ax = gca; hold all;
% figure; off_ax = gca; hold all

c = colormap('Lines');

%do hulls
type_hull = cell(num_types,1);
hull_edges = cell(num_types,1);
for k = 1:num_types   
    for c = C.type.(types{k})
        cell_dat = cell_data(c);
        h = [];
        [h(:,1), h(:,2)] = poly2cw(cell_dat.hull_2d(:,1), cell_dat.hull_2d(:,2));
        if isempty(type_hull{k})
            type_hull{k} = h;
        else
            temp_hull = [];
            [temp_hull(:,1), temp_hull(:,2)] = polybool('union', type_hull{k}(:,1), type_hull{k}(:,2), h(:,1), h(:,2));
            type_hull{k} = temp_hull;
        end
    end
    
    hull_edges{k} = zeros(size(type_hull{k}));
    last_beginning = 1;
    current_node = 1;
    remove_this_line = false(size(hull_edges{k},1),1);
    for n = 1:size(type_hull{k},1)
        if isnan(type_hull{k}(n,1))
            hull_edges{k}(n-1,2) = last_beginning;
            remove_this_line(n) = true;
            last_beginning = n+1;
        else
            hull_edges{k}(n,:) = [current_node current_node+1];
        end
        current_node = current_node+1;
    end
    hull_edges{k}(end,2) = last_beginning;
    hull_edges{k}(remove_this_line,:) = [];
end

sac_content_matrix = zeros(num_types, num_sacs, num_bins);
sac_contact_matrix = zeros(num_types, num_sacs, num_bins);

c_types = zeros(1000,1);
c_nums = zeros(1000,1);
c = 0;
for k = 1:num_types
    cns = C.type.(types{k});
    c_nums(c+(1:length(cns))) = cns;
    c_types(c+(1:length(cns))) = k;    
    c = c+length(cns);
end
c_nums = c_nums(1:c);
c_types = c_types(1:c);

cell_num2type = sparse(c_nums,ones(length(c_nums),1),c_types,100000,1);


for s = 1:num_sacs
    sac_dat = cell_data(sac_nums(s));
    sac_mid = sac_dat.get_midpoint(true);
    
    p = sac_dat.get_surface;
    depth = C.f(p(:,1));
    p(depth < 15,:) = [];
    
    
    d = sqrt((p(:,2)-sac_mid(2)).^2 + (p(:,3)-sac_mid(3)).^2);
    for k = 1:num_types
        %is_in_hull = inpolygon(p(:,2),p(:,3), type_hull{k}(:,1),type_hull{k}(:,2));
        is_in_hull = inpoly(p(:,2:3), type_hull{k}, hull_edges{k});
        sub_d = d(is_in_hull)/1000;
    
        for n = 1:num_bins
            sac_content_matrix(k,s,n) = sum(sub_d>=bins(n)-bin_size/2 & sub_d<bins(n)+bin_size/2);
        end
    end
    
    
    conts = double(sac_dat.contacts);
    depth = C.f(conts(3,:));
    conts = conts(:,depth>15);
    
    cont_type = full(cell_num2type(conts(1,:)))';
    d = sqrt((conts(4,:) - sac_mid(2)).^2 + (conts(5,:) - sac_mid(3)).^2)/1000;
    bin_num = ceil((d-bins(1)+bin_size/2)/bin_size);
    
    for k = 1:num_types
        for n = 1:num_bins
            sac_contact_matrix(k,s,n) = sum(conts(2, bin_num==n & cont_type==k));
        end
    end
    
end

M = zeros(numel(sac_contact_matrix),5);
n = 1;
for s = 1:num_sacs
    for k = 1:num_types
        for b = 1:num_bins
            M(n,1) = sac_nums(s);
            M(n,2) = k;
            M(n,3) = b;
            M(n,4) = sac_contact_matrix(k,s,b);
            M(n,5) = sac_content_matrix(k,s,b);
            n=n+1;
        end
    end
end

csvwrite('supp_fig_tot.txt', M);

sac_contact_matrix = sac_contact_matrix*291.5/1000/1000;
sac_content_matrix = sac_content_matrix*1.4*291.5/1000/1000;



x_ax_vals = bins;

figure; hold on
for k = 1:5;
    plot(x_ax_vals, squeeze(sum(sac_contact_matrix(k,:,:),2))', 'LineWidth', 2, 'Color', C.colormap(k,:));
end
prep_figure(gcf, gca, ...
    'xlabel', 'Distance from soma', ...
    'ylabel', 'Total Contact (um^2)' ...
);

figure; hold on
for k = 1:5;
    plot(x_ax_vals, squeeze(sum(sac_content_matrix(k,:,:),2))', 'LineWidth', 2, 'Color', C.colormap(k,:));
end
prep_figure(gcf, gca, ...
    'xlabel', 'Distance from soma', ...
    'ylabel', 'SAC surface area' ...
);

sac_dens_matrix = squeeze(sum(sac_contact_matrix,2)./sum(sac_content_matrix,2))'*100;
figure; hold on

for k = 1:5
    plot(x_ax_vals, sac_dens_matrix(:,k), 'LineWidth', 2, 'Color', C.colormap(k,:));
end

prep_figure(gcf, gca, ...
    'xlabel', 'Distance from soma', ...
    'ylabel', 'SAC coverage (percent)' ...
    );

        
    
% legend({'type 1', 'type 3a'})
% legend({'BC1', 'BC2', 'BC3a', 'BC3b', 'BC4'})
% ylabel('BC contact per SAC area (um^2/um^2)', 'FontSize', 20)
% xlabel('Distance from soma (microns)', 'FontSize', 20)
% 
% set(gca, 'YTick', [0:.5:10]);
% 
% set(gcf, 'Position', [0 0 640 640]);
% set(gca, 'FontSize', 20);
% 
% set(gca, 'Position', [.1 .1 .8 .8]);




