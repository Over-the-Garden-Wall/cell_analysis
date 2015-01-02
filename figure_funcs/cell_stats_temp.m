function M = cell_stats(cell_nums)
    C = get_constants;

    for n = 1:length(cell_nums)
        percentile_depths = [.1 .2 .3 .4 .5 .6 .7 .8 .9];
        
        c = cell_nums(n);
        c_d = cell_data(c);
        
        p = c_d.get_surface;
        d = C.f(p(:,1));
        
        total_ps = length(d);
        
        pds = d(round(percentile_depths*total_ps));
        
        s = c_d.stratification;
        
        
        trunk_thickness =  sum(d>10 & d< 30);
        
        [mode_val mode_strat] = max(s);
        
        trunk_end = find(s(1:mode_strat) < .02, 1, 'last');
        trunk_end = C.strat_x(trunk_end);
        
        mode_strat = C.strat_x(mode_strat);
        
        d(d<trunk_end) = [];
        mse = mean((d - mode_strat).^2);
        
        M(n,:) = [c, total_ps, trunk_thickness, c_d.hull_area/10^6, trunk_end, mode_strat, mse, pds'];
        
    end
    
end