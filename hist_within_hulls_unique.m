
C = get_constants;

types = {'t1', 't2', 't3a', 't3b', 't4'};

ref_type = 'sure_off_sac';


bin_size = 15;
bins = bin_size/2:bin_size:bin_size*8.5;


num_types = length(types);
num_bins = length(bins);

hulls = cell(num_types,1);
cell_mids = cell(num_types,1);

for k = 1:num_types
    hulls{k} = [];
    for c = 1:length(C.type.(types{k}))
        cell_dat = cell_data(C.type.(types{k})(c));
        temp = [];
        [temp(:,1), temp(:,2)] = poly2cw(cell_dat.hull_2d(:,1), cell_dat.hull_2d(:,2));
        if c > 1
            temp2 = [];
            [temp2(:,1), temp2(:,2)] = polybool('union', temp(:,1), temp(:,2), hulls{k}(:,1), hulls{k}(:,2));
        else
            temp2 = temp;
        end
        hulls{k} = temp2;        
    end
end

bin_hist = zeros(num_bins,num_types);
for s = C.type.(ref_type);
    cell_dat = cell_data(s);
    soma = cell_dat.get_midpoint(true);
    p = cell_dat.get_surface;
    
    d = C.f(p(:,1));
    p = p(d >= 15 & d<= 70,2:3);
    d = sqrt((p(:,1) - soma(2)).^2 + (p(:,2) - soma(3)).^2)/1000;
    bin_num = ceil(d/bin_size);
        
    for k = 1:num_types
        is_in_hull = rough_2d_inpolygon(p, hulls{k}, [200 200]);
        for n = 1:num_bins
            bin_hist(n,k) = bin_hist(n,k)+sum(is_in_hull(bin_num==n));            
        end
        
    end
    
    
end

figure; 
hold on; 
for k = 1:5; 
    h(k) = plot(bins,bin_hist(:,k)*410/1000/1000, 'LineWidth', 2, 'Color', C.colormap(k,:)); 
end

prep_figure(gcf, gca, 'xlabel', 'Distance from soma (um)', 'ylabel', 'Surface area (um^2)')

% plot(bins,bin_hist*410/1000/1000/length(C.type.(ref_type)), 'LineWidth', 2, 'Color', C.colormap(6,:));
% figure; plot(bins,bin_hist, 'LineWidth', 2);