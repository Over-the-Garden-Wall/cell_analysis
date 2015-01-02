function prep_figure(varargin)

    ip = inputParser;
    ip.addOptional('figure_handle',gcf, @(x) ishandle(x));
    ip.addOptional('axis_handle',gca, @(x) ishandle(x));
    ip.addParamValue('xlabel', [], @(x) ischar(x));
    ip.addParamValue('ylabel', [], @(x) ischar(x));
    ip.addParamValue('title', [], @(x) ischar(x));
    ip.addParamValue('legend', [], @(x) iscell(x));
    ip.addParamValue('legend_handles', [], @(x) ishandle(x));    
    ip.addParamValue('fontsize', 12/.45, @(x) isnumeric(x));
    ip.addParamValue('YTick', [], @(x) isnumeric(x));
    ip.addParamValue('XTick', [], @(x) isnumeric(x));    
    ip.addParamValue('size', [640 640], @(x) isnumeric(x));
    
    ip.parse(varargin{:});
    s = ip.Results;
    
    if ~isempty(s.xlabel)
        xlabel(s.axis_handle, s.xlabel, 'FontSize', s.fontsize);
    end
    if ~isempty(s.xlabel)
        ylabel(s.axis_handle, s.ylabel, 'FontSize', s.fontsize);
    end
    if ~isempty(s.XTick)
        set(s.axis_handle, 'XTick', s.XTick);
    end
    if ~isempty(s.YTick)
        set(s.axis_handle, 'YTick', s.YTick);
    end
    if ~isempty(s.title)
        title(s.title, 'fontSize', 20);
    end
    if ~isempty(s.legend)
        if isempty(s.legend_handles)            
            legend(s.axis_handle, s.legend, 'FontSize', s.fontsize)
        else
            legend(s.axis_handle, s.legend_handles, s.legend, 'FontSize', s.fontsize)
        end
        legend(s.axis_handle, 'boxoff');
    end
    
    
    set(s.figure_handle, 'Position', [0 0 s.size]);
    
    
    set(s.axis_handle, 'FontSize', 20);

    set(s.axis_handle, 'Position', [.1 .1 .8 .8]);

    set(s.axis_handle, 'box', 'off');

end