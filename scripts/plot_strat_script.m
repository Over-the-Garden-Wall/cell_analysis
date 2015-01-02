function plot_strat_script(cells, target_cell, varargin)
% 
% target_cell is the cell with which contacts will be annotated. 
% It defaults to 10010 if empty or plot_strat_script is called with 1 argument.
% 
% the cells variable can either be a vector of cell ids, and string, or a cell array.
% 
% if cells is a vector of cell ids, it will plot the stratification of those cells.
% 
% if cells is a string, it will check if the string is 'j', 'off_sac', 't1', 
% 't2', 't3a', 't3b', or 't4'. If it is none of those, there will be an error.
% Otherwise, it will plot the cells that belong to that group. The list of 
% cell types can be seen in get_constants.m, and correspond to sebastian's recent designation
% 
% if cells is a cell array, then each element must be a vector of cell ids or string as above. 
% All appropriate cells will be plotted, and colored according to which element they belong to.
% 
% Example 1: plot_strat_script('t1') will plot the type 1 bipolar cells
% Example 2: plot_strat_script([60008 60012 60019 60026]) will plot these 4 cells
% Example 3: plot_strat_script({[60008 60012 60019 60026], 't1'}) will plot these 4 
% cells in one color, and all the type 1 bipolar cells in a different color




    if ~exist('target_cell','var') || isempty(target_cell)
        target_cell = 10010;
    end

    C = get_constants;

    if ~iscell(cells)
        temp = cells;
        clear cells;
        cells{1} = temp;
    end
    
    cell_nums = [];
    cell_type = [];
    for k = 1:length(cells)
        if isstr(cells{k})
            fns = fieldnames(C.type);
            fn = 0;
            for f = 1:length(fns);
                if strcmp(fns{f},cells{k})
                    fn = f;
                end
            end
            cell_nums = [cell_nums C.type.(fns{fn})];
            cell_type = [cell_type k*ones(1,length(C.type.(fns{fn})))];
        else
            cell_nums = [cell_nums cells{k}];
            cell_type = [cell_type k*ones(1,length(cells{k}))];            
        end
    end

    conn_data = load(C.conn_loc);
    fns = fieldnames(conn_data);
    
    if isempty(varargin)
        plot_strat_w_contacts(cell_nums, target_cell, conn_data.(fns{1}), cell_type);
    else
        plot_strat_w_contacts(cell_nums, target_cell, conn_data.(fns{1}), varargin{:});
    end
end