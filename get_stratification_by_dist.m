function [quartile_data bins full_x full_data] = get_stratification_by_dist(cell_num, varargin)
    
    C = get_constants;
    
    num_cells = length(cell_num);
    
    
    
    p = inputParser;    
    p.addRequired('cell_num', @isnumeric);
    p.addOptional('dist_axis', [], @(x) isempty(x) || ((isnumeric(x) && size(x,2) == 2) && (size(x,1)==1 || size(x,1)==num_cells)));
    p.addOptional('num_bins', 10, @isnumeric);
    p.addOptional('quartile_points', [.00001 .1 .25 .5 .75 .9 1], @isnumeric);
    p.addOptional('use_soma', [], @(x) isempty(x) || islogical(x))
    p.addOptional('rescale_individually', false, @islogical)
    p.addOptional('x_step', 1, @isnumeric)
        
    
    p.parse(cell_num, varargin{:});
    
    s = p.Results;
    
    default_axis = zeros(num_cells,2);
    for n = 1:num_cells
        cell_dat{n} = cell_data(cell_num(n));
        default_axis(n,:) = cell_dat{n}.dist_axis;
    end
    
    if isempty(s.dist_axis)
        s.dist_axis = default_axis;
        s.dist_axis = zeros(size(default_axis));
    end
    
    
    if size(s.dist_axis,1) == 1;
        s.dist_axis = ones(num_cells,1)*s.dist_axis;
    end
    
    
    
    dist = [];
    depth = [];
    for n = 1:num_cells
        
        cell_mid = cell_dat{n}.get_midpoint(s.use_soma);

        cell_p = cell_dat{n}.get_surface;

        for d = 2:3
            cell_p(:,d) = cell_p(:,d)-cell_mid(d);
        end

        if ~all(s.dist_axis == 0)
            my_dist = cell_p(:,2) * s.dist_axis(n,1) + cell_p(:,3) * s.dist_axis(n,2);            
        else
            my_dist = sqrt(cell_p(:,2).^2 + cell_p(:,3).^2);        
        end
        
        if s.rescale_individually
            my_dist = my_dist/max(my_dist)*100;
        else
            my_dist = my_dist/1000;
        end
        
        
        dist = [dist; my_dist];
%         warning('hack in place');
%         disp(max(cell_p(:,2) * s.dist_axis(n,1) + cell_p(:,3) * s.dist_axis(n,2)));
        
        %warning('HACKKKKKK!!!!')
%         depth = [depth; cell_p(:,2) * -s.dist_axis(n,2) + cell_p(:,3) * s.dist_axis(n,1)];

        depth = [depth; C.f(cell_p(:,1))];
        
        
    end
    
    
    dist_min = min(dist);
    dist_max = max(dist);
    
    if length(s.num_bins) == 1
        bin_size = (dist_max-dist_min)/s.num_bins;
        bins = (1:s.num_bins)*bin_size - bin_size/2 + dist_min;
    else
        bins = s.num_bins;
        bin_size = bins(2)-bins(1);
        if any(bins(2:end)-bins(1:end-1) ~= bin_size)
            error('bins are of unequal size');
        end
        s.num_bins = length(bins);
    end
        
    num_quarts = length(s.quartile_points);
    quartile_data = zeros(num_quarts, s.num_bins);
    
    min_depth = floor(min(depth));
    full_x = min_depth:s.x_step:ceil(max(depth));
    full_data = zeros(length(full_x), s.num_bins);
        
    for n = 1:s.num_bins
        is_in_bin = (dist >= bins(n) - bin_size/2) & (dist <= bins(n) + bin_size/2);
        if any(is_in_bin)
            full_data(:,n) = hist(depth(is_in_bin), full_x);

            data_cumsum = cumsum(full_data(:,n));
            data_cumsum = data_cumsum/data_cumsum(end);

            for k = 1:num_quarts            
                quartile_data(k,n) = full_x(find(data_cumsum >= s.quartile_points(k),1,'first'));
            end
        else
            quartile_data(:,n) = NaN;
        end
    end
end
        