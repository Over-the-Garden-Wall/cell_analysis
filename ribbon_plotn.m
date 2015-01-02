function [axes_hand line_hand] = ribbon_plot2(x_cell, y_cell, ribbon_size, varargin)

    p = inputParser;    
    p.addRequired('x_cell', @iscell);
    p.addRequired('y_cell', @iscell);
    p.addRequired('ribbon_size', @iscell);
    
    p.addParamValue('x_spacing', .01, @isnumeric);
    p.addParamValue('title', [], @ischar);
    p.addParamValue('has_legend', true, @islogical);
    p.addParamValue('colors', colormap('Lines'), @(x)size(x,2)==3);
    p.addParamValue('alpha', .2, @(x) isnumeric(x) && x>=0 && x<=1);
    p.addParamValue('edgeVal', 'zero', @(x) strcmp(x,'zero') || strcmp(x,'nearest'));
    
    
    
    p.parse(x_cell, y_cell, ribbon_size, varargin{:});
    
    s = p.Results;
    
%     axes_hand = gca;
    
    
    num_lines = length(x_cell);
    
    all_x = [];
    for k = 1:num_lines
        all_x = [all_x; x_cell{k}(:)];
    end
    all_x = unique(all_x);
    
    
    num_points = length(all_x);            
    
    all_y = zeros(num_points, num_lines);
    all_rib = zeros(num_points, num_lines);
    rib_ends = zeros(num_points, num_lines, 2);
    
    interp_with_zero = strcmp(s.edgeVal,'zero');
    
    for k = 1:num_lines
        all_y(:,k) = interp1(x_cell{k}, y_cell{k}, all_x);   
        all_rib(:,k) = interp1(x_cell{k}, ribbon_size{k}, all_x);   
        
        start_point = find(~isnan(all_y(:,k)),1,'first');
        end_point = find(~isnan(all_y(:,k)),1,'last');
        
        if interp_with_zero
            all_y([1:start_point-1, end_point+1:end],k) = 0;
            all_rib([1:start_point-1, end_point+1:end],k) = 0;            
        else
            all_y(1:start_point-1,k) = all_y(start_point,k);
            all_y(end_point+1:end,k) = all_y(end_point,k);
            all_rib(1:start_point-1,k) = all_rib(start_point,k);
            all_rib(end_point+1:end,k) = all_rib(end_point,k);
        end
    
        rib_ends(:,:,1) = all_y - all_rib;
        rib_ends(:,:,2) = all_y + all_rib;
        
    end
    
    intersect_P = [];
    for k = 1:size(rib_ends(:,:),2)
        for l = k+1:size(rib_ends(:,:),2)
            intersect_P = [intersect_P, InterX([all_x'; rib_ends(:,k)'], [all_x'; rib_ends(:,l)'])];
        end
    end
    
    all_x = [all_x; intersect_P(1,:)'];
    
    all_x = unique(all_x);
    
    
    num_points = length(all_x);            
    
    all_y = zeros(num_points, num_lines);
    all_rib = zeros(num_points, num_lines);
    rib_ends = zeros(num_points, num_lines, 2);
    
    interp_with_zero = strcmp(s.edgeVal,'zero');
    
    for k = 1:num_lines
        all_y(:,k) = interp1(x_cell{k}, y_cell{k}, all_x);   
        all_rib(:,k) = interp1(x_cell{k}, ribbon_size{k}, all_x);   
        
        start_point = find(~isnan(all_y(:,k)),1,'first');
        end_point = find(~isnan(all_y(:,k)),1,'last');
        
        if interp_with_zero
            all_y([1:start_point-1, end_point+1:end],k) = 0;
            all_rib([1:start_point-1, end_point+1:end],k) = 0;            
        else
            all_y(1:start_point-1,k) = all_y(start_point,k);
            all_y(end_point+1:end,k) = all_y(end_point,k);
            all_rib(1:start_point-1,k) = all_rib(start_point,k);
            all_rib(end_point+1:end,k) = all_rib(end_point,k);
        end
    
        rib_ends(:,:,1) = all_y - all_rib;
        rib_ends(:,:,2) = all_y + all_rib;
    
        rib_ends(isnan(rib_ends)) = 0;
        all_y(isnan(all_y)) = 0;
        
        
    end
    
    
    %let's do this the stupid way, (2^n)*n  areas
    num_areas = 2^num_lines * num_lines;
    region_areas = zeros(num_points, num_areas);
    
    region_ids = false(num_lines,num_areas);
    for k = 1:num_lines
        region_ids(k,mod(floor((0:(num_areas-1))/2^(k-1)),2)==1) = true;
    end
        
    
    
    
    for x = 1:num_points
        
        old_region = 1;
        region = false(num_lines,1);
        old_y = 0;
        
        relevant_y = unique(rib_ends(x,:));
        
        for y = relevant_y
            for k = 1:num_lines
                if y == rib_ends(x,k,2)
                    region(k) = false;
                elseif y == rib_ends(x,k,1)
                    region(k) = true;
                end
            end
            
            
            region_areas(x, old_region) = y-old_y;
            old_y = y;
            
            k = old_region;
            while 1
                k = k+1;
                if k > 1+num_areas
                    error('sadface');
                end
                if k == 1+num_areas || all(region_ids(:,k)==region)
                    break
                end
            end
            old_region = k;
            
        end
                            
    end
    

    region_colors = zeros(num_areas,3);
 
    
    for k = 1:num_areas
        num_trues = sum(region_ids(:,k));
        if num_trues == 0
            region_colors(k,:) = [1 1 1];
        else
            
            for n = 1:num_lines
                if region_ids(n,k)
                    region_colors(k,:) = region_colors(k,:) + s.colors(k,:)/num_trues;
                end
            end
        end
    end                        
            
    region_colors = region_colors*s.alpha+(1-s.alpha);
    
%     set(gcf,'DefaultAxesColorOrder',[region_colors; s.colors]);
    
%     set(gcf,'DefaultAxesColorOrder',[region_colors; s.colors]);
    axes_hand = gca;
    hold(gca, 'all');
    
    area_hands = area(all_x, region_areas,'lineStyle', 'none');
    
    for k = 1:num_lines
        line_hand(k) = plot(x_cell{k}, y_cell{k});
    end
    colormap(region_colors);
    
%     set(area_hands(1), 'Visible', 'off');
%     set(area_hands(5), 'Visible', 'off');
    
    
    title(s.title);
    
end
    
    
            
        
        
        
