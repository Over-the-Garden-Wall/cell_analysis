function [axes_hand line_hand] = ribbon_plot2(x_cell, y_cell, ribbon_size, varargin)

    p = inputParser;    
    p.addRequired('x_cell', @iscell);
    p.addRequired('y_cell', @iscell);
    p.addRequired('ribbon_size', @iscell);
    
    p.addParamValue('x_spacing', .01, @isnumeric);
    p.addParamValue('title', [], @ischar);
    p.addParamValue('has_legend', true, @islogical);
    p.addParamValue('colors', [1 0 0; 0 0 1], @(x)size(x,2)==3);
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
    for k = 1:4
        for l = k+1:4
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
    
    %assuming r and b are the colors, area coding is: w r b p w b r
    region_areas = zeros(num_points,8);
    
    for x = 1:num_points
        region_areas(x,1) = min(rib_ends(x,:));
        if rib_ends(x,1,1) < rib_ends(x,2,1)
            lower_num = 1;
        else
            lower_num = 2;
        end
        higher_num = 3-lower_num;
        
        if rib_ends(x,lower_num,2) < rib_ends(x,higher_num,1)
            has_intersection = false;
        else
            has_intersection = true;
        end
        
        if has_intersection
            region_areas(x,1+lower_num) = rib_ends(x,higher_num,1) - rib_ends(x,lower_num,1);
            region_areas(x,4) = rib_ends(x, lower_num, 2) - rib_ends(x,higher_num,1);
            region_areas(x,5+higher_num) = rib_ends(x,higher_num,2) - rib_ends(x,lower_num,2);
        else
            region_areas(x,1+lower_num) = rib_ends(x,lower_num,2) - rib_ends(x,lower_num,1);
            region_areas(x,5) = rib_ends(x, higher_num, 1) - rib_ends(x,lower_num,2);
            region_areas(x,5+higher_num) = rib_ends(x,higher_num,2) - rib_ends(x,higher_num,1);
        end
        
        rib_ends(isnan(rib_ends)) = 0;
        all_y(isnan(all_y)) = 0;
        
    end
    
    region_colors = [1 1 1; s.colors(2,:); s.colors(1,:); s.colors(1,:)/2 + s.colors(2,:)/2; ...
        1 1 1; s.colors(2,:); s.colors(1,:)];
%     region_colors = [ s.colors(2,:); s.colors(1,:); s.colors(1,:)/2 + s.colors(2,:)/2; ...
%          s.colors(2,:); s.colors(1,:)];
    
    region_colors = region_colors/4+3/4;
    
%     set(gcf,'DefaultAxesColorOrder',[region_colors; s.colors]);
    
%     set(gcf,'DefaultAxesColorOrder',[region_colors; s.colors]);
    axes_hand = gca;
    hold(gca, 'all');
    
    area_hands = area(all_x, region_areas,'lineStyle', 'none');
    
    line_hand = plot(all_x, all_y(:,1), all_x, all_y(:,2));
    colormap(region_colors);
    
%     set(area_hands(1), 'Visible', 'off');
%     set(area_hands(5), 'Visible', 'off');
    
    
    title(s.title);
    
end
    
    
            
        
        
        
