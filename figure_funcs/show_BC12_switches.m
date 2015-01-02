function show_BC12_switches

    C = get_constants;   

    BCs = [C.type.t1, C.type.t2];
    actual_type = [ones(1,length(C.type.t1)) 2*ones(1,length(C.type.t2))];
    
    
    
    num_BCs = length(BCs);
    
    %need 10th, 25th, 75, percentiles and axonal volume
    split_vals = zeros(num_BCs,4);
    
    for n = 1:num_BCs
        cell_dat = cell_data(BCs(n));
        p = cell_dat.get_surface;
        d = C.f(p(:,1));
        
        d(d<0) = [];
        
        d = sort(d);
        num_points = length(d);
        split_vals(n,1:3) = d(ceil([.1 .25 .75]*num_points));
        
        split_vals(n,4) = cell_dat.V;
    end
    
    split = split_vals(:,3) - split_vals(:,2);
    pred_type = 1+(split' >= 10.3);
    h = zeros(3,1);
    figure; hold on
    mismatch = false(1,size(split_vals,1));
    for k = 1:2
        h(k) = scatter(split_vals(actual_type==k,2), split_vals(actual_type==k,3), '*', 'MarkerEdgeColor', C.colormap(k,:));
        mismatch = mismatch | (actual_type==k & pred_type ~=k);
        
    end
    
    h(3) = scatter(split_vals(mismatch,2), split_vals(mismatch,3), 'o', 'MarkerEdgeColor', [0 0 0], 'sizeData', 72);
    
    prep_figure(gcf,gca, 'xlabel', 'first quartile depth (% IPL)', 'ylabel', 'third quartile depth (% IPL)', 'legend', {'BC1', 'BC2', 'Swapped'});
    set(gca, 'XTick', 0:5:30);
    set(gca, 'YTick', 0:5:40);
    
end
    
    