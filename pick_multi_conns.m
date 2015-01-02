function selected_R = pick_multi_conns(R, cells, varargin)
    
    p = inputParser;    
    p.addRequired('R', @isnumeric);
    p.addRequired('cells', @(x)iscell(x) && length(x)==2 && isnumeric(x{1}));    
    p.addOptional('cell_output', true, @islogical);
    
    p.parse(R, cells, varargin{:});
    s = p.Results;
    
    if s.cell_output
        selected_R = cell(length(s.cells{1}), length(s.cells{2}));
        for k = 1:length(s.cells{1})
            selected_R(k,:) = pick_conns(s.R, s.cells{1}(k), s.cells{2}, s.cell_output);
        end
    else
        selected_R = [];
        for k = 1:length(cells{1})
            selected_R = [selected_R pick_conns(s.R, s.cells{1}(k), s.cells{2}, s.cell_output)];
        end
    end
end
    