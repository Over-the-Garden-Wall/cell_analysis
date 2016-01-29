function pr_pred = peters_rule_prediction(ref_cells, comp_cells, depth_range)

    if ~exist('depth_range', 'var') || isempty(depth_range)
        depth_range = [0 100];
    end

    pr_pred = zeros(length(comp_cells),1);
    ref_strat = get_gridspace_strat(ref_cells, depth_range, false);
    
    for n = 1:length(comp_cells)
        comp_strat = get_gridspace_strat(comp_cells{n}, depth_range, true);
        pr_pred(n) = sum(ref_strat.*comp_strat)/sum(ref_strat);
    end
    
end
   
function strat_per_cubicnm = get_gridspace_strat(cell_nums, depth_range, use_volume)

    C = get_constants;

    total_strat = zeros(depth_range(2)-depth_range(1)+1,1);
    total_hull = [];
    
    for c = cell_nums;
        c_d = cell_data(c);
        s = c_d.stratification;
        s = s(1+depth_range(1)-C.strat_x(1):end);
        s(end+1:size(total_strat,1)) = 0;
        if use_volume
            total_strat = total_strat + s(1:size(total_strat,1))*c_d.V;
        else
            total_strat = total_strat + s(1:size(total_strat,1))*c_d.SA;
        end
        
        h = cell(1,2);
        [h{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
            
        if isempty(total_hull)
            total_hull = h;
        else
            [total_hull{:}] = polybool('union', h{1}, h{2}, total_hull{1}, total_hull{2});            
        end
            
    end
    total_hull_area = poly_area([total_hull{1}, total_hull{2}]);
       
    perc_change_for_nm = C.f([1; 2]);
    perc_height = abs(1/(perc_change_for_nm(2) - perc_change_for_nm(1)));
    
    strat_per_cubicnm = total_strat * prod(C.res) / total_hull_area / perc_height;
end
    