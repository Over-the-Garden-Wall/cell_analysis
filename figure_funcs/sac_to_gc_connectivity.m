function gc_nums = sac_to_gc_connectivity
    
    C = get_constants;
    gc_nums = C.type.ganglion;
    
    sac_nums = {C.type.on_sac C.type.off_sac};
    l_lims = [22 42; 52 72];
    
    sac_conns = cell(2,1);                
    
    for l = 1:2
        sac_conns{l} = cell(length(sac_nums{l}),1);
        for sacn = 1:length(sac_nums{l})
            s_conts{l}{sacn} = detect_vericose_contacts(sac_nums{l}(sacn), 500, 200, 50000);
        end
    end
    
    %greedy cluster   
    gc_strats = zeros(100,length(gc_nums));
    for n = 1:length(gc_nums)
        try
            c_d = cell_data(gc_nums(n));
            s = c_d.stratification;
            s = s(1-C.strat_x(1):end);        
            gc_strats(1:min(length(s),100), n) = s(1:min(100,length(s)));        
        catch ME
            disp(ME.message)
        end
    end    
    for n = 1:length(gc_nums)-1
        strat_sim = -ones(length(gc_nums),1);
        for k = n+1:length(gc_nums)
            strat_sim(k) = gc_strats(:,n)'*gc_strats(:,k);
        end
        [~, next_gc] = max(strat_sim);
        gc_nums([n+1, next_gc]) = gc_nums([next_gc, n+1]);
        gc_strats(:, [n+1, next_gc]) = gc_strats(:, [next_gc, n+1]);
    end
    
    
    
    gc_area = zeros(length(gc_nums),2);
    gc_sac_conn = zeros(length(gc_nums),2);
    for l = 1:2
        for gcn = 1:length(gc_nums)
            gc = gc_nums(gcn);
            try    
                for sacn = 1:length(sac_nums{l})
                    gc_sac_conn(gcn,l) = gc_sac_conn(gcn,l) + sum(s_conts{l}{sacn}(2, s_conts{l}{sacn}(1,:) == gc));
                end


                gc_d = cell_data(gc);
                p = gc_d.get_surface;
                d = C.f(p(:,1));
                p = p(d >= l_lims(l, 1) & d <= l_lims(l, 2),:);
            
                p_hull = convhull(p(:,2), p(:,3));
                p = p(p_hull, 2:3);
                [p(:,1), p(:,2)] = poly2cw(p(:,1), p(:,2));
            
                gc_area(gcn, l) = polyarea(p(:,1), p(:,2));
            catch ME
                disp(ME.message)
            end
        end
    end
            
    figure; hold all
    for l = 1:2
        scatter((1:length(gc_nums))', gc_sac_conn(:,l)./gc_area(:,l));
    end
    set(gca, 'YScale', 'log')
    
end
            
            
            
    
    
    