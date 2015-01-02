function [density is_valid] = get_density_by_depth(cell_num, force_recalc)

    if ~exist('force_recalc','var') || isempty(force_recalc)
        force_recalc = false;
    end

    C = get_constants;
    
    strat_fn = [C.strat_dir '/cell_' num2str(cell_num) '_strat.mat'];
    if exist(strat_fn,'file') && ~force_recalc
        load(strat_fn);
        is_valid = true;
    else    
    
        x = C.strat_x;
        
        fn = [C.point_dir '/cell_' num2str(cell_num) '_surface.mat'];
        if exist(fn,'file')
            data = load(fn);
            fields = fieldnames(data);
            p = C.f(data.(fields{1})(:,1));
            density = hist(p,x);            
            is_valid = true;
            save(strat_fn, 'density');
        else
            warning(['data not found for ' fn]);
            density = zeros(size(x));        
            is_valid = false;
        end

    end
    
end