function [line_handles area_handles fig_handle] = plot_with_ste_area(plot_x, plot_y, ribbon_size, fig_handle)

    if ~iscell(plot_x)
        temp = plot_x;
        clear plot_x
        plot_x{1} = temp;
    end
    if ~iscell(plot_y)
        temp = plot_y;
        clear plot_y
        plot_y{1} = temp;
    end
    if ~iscell(ribbon_size)
        temp = ribbon_size;
        clear ribbon_size
        ribbon_size{1} = temp;
    end

    if ~exist('fig_handle', 'var') || isempty(fig_handle)
        fig_handle = gcf;
        axes_h = gca;
    else
        axes_h = get(fig_handle, 'CurrentAxes');
        if isempty(axes_h)
            axes_h = gca;
        end        
    end
    
    hold(axes_h,'all');
    
    num_areas = length(plot_x);
    
    area_handles = zeros(num_areas,2);
    line_handles = zeros(num_areas,1);
    
    
%     x = -9:80;

    for n = 1:num_areas        
        if size(plot_x{n},1) ~= 1;
            plot_x{n} = plot_x{n}';
        end
        if size(plot_y{n},1) ~= 1;
            plot_y{n} = plot_y{n}';
        end
        if size(ribbon_size{n},1) ~= 1;
            ribbon_size{n} = ribbon_size{n}';
        end                
        
        
        area_handles(n,:) = area(axes_h, plot_x{1}, [plot_y{n}-ribbon_size{n}/2; ribbon_size{n}]', 'lineStyle', 'none');
        set(area_handles(1), 'Visible', 'off');
        line_handles(n) = plot(axes_h, plot_x{n}, plot_y{n});
        line_col = get(line_handles(n), 'Color');
        set(area_handles(2), 'FaceColor', line_col*.5 + [1 1 1]*.5);
    
    end
    
    
%     a = 1;
end