load('sac2j1conns.mat');

all_dist = sac2j1conns.dist;
all_cont = sac2j1conns.area;
cell_nums = sac2j1conns.cell_id;

[all_dist sort_ord] = sort(all_dist/1000);
all_cont = all_cont(sort_ord);
cell_nums = cell_nums(sort_ord);

    num_bins = 5;
    
    
    x = zeros(num_bins,1);
    y = zeros(num_bins,1);
    ste_y = zeros(num_bins,1);
    
    n_per_bin = length(all_dist)/num_bins;
    
    for n = 1:num_bins
        valid_range = floor((n-1)*n_per_bin)+1:floor(n*n_per_bin);
        
        x(n) = mean(all_dist(valid_range));
        y(n) = mean(all_cont(valid_range));
        ste_y(n) = std(all_cont(valid_range)) / sqrt(length(valid_range)-1);
        
    end
    
    %
    %         x = 1:C.dist_bin:max(all_dist);
    %         y = zeros(size(x));
    %         ste_y = zeros(size(x));
    %
    %         for k = 1:length(x)
    %
    %             valid_dist = all_dist>=x(k) & all_dist<x(k)+C.dist_bin;
    %             if any(valid_dist)
    %                 y(k) = mean(all_cont(valid_dist));
    %                 ste_y(k) = std(all_cont(valid_dist)) / sqrt(sum(valid_dist)-1);
    %             end
    %         end
    %         x = x + C.dist_bin/2;
    %
    ste_y(isnan(ste_y)) = 0;
    
% end

figure;
[~, line_hands] = ribbon_plotn({x}, {y}, {ste_y});


title( 'Connectivity of Off-SACs with J-Cell')
ylabel( 'Contact Area per Synapse')
xlabel( 'Distance along directionally selective axis (microns from soma)');
% legend(line_hands, {'Type 1/2 Bipolar Cells', 'Type 3/4 Bipolar Cells'});

gca; hold(gca,'all');
% for k = 1:2
    scatter(all_dist, all_cont, '*');
% end
%         for n = 1:length(all_dist)
%
% %             text(all_dist(n), all_cont(n) + max(all_cont)/30, num2str(cell_nums(n)));
% %             SA = get_size_stats(cell_nums(n));
% %             all_cont(n) = all_cont(n)/SA;
%         end

% legend({'Type 1/2 to J contacts', 'Type 3/4 to J contacts'});

figure;
[hist_y x_out] = hist(all_dist,10);
plot(x_out, hist_y);

% title( 'contact area by distance from soma scatter')
% ylabel( 'Contact Area per Synapse')
% xlabel( 'distance along soma-distal axis (microns)');
