function show_BC12_switches

    C = get_constants;   

    BCs = [C.type.t3a, C.type.t3b, C.type.t4];
    actual_type = [ones(1,length(C.type.t3a)) 2*ones(1,length(C.type.t3b)) 3*ones(1,length(C.type.t4))];
    
    
    
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
    
    
    is_4 = (split_vals(:,1) < 24);
    is_3 = (split_vals(:,1) >= 24);

    is_3a = (split_vals(:,4) >= 9.09*10^5) & is_3;
    is_3b = (split_vals(:,4) < 9.09*10^5) & is_3;
    
    
    
%     is_3a = (split_vals(:,4) >= 9.5*10^5);
%     is_3b4 = (split_vals(:,4) < 9.5*10^5);
%     
%     is_4 = (split_vals(:,1) < 22) & is_3b4;
%     is_3b = (split_vals(:,1) >= 22) & is_3b4;
%     
    pred_type = is_3a' + is_3b'*2 + is_4'*3;
    
    
    h = zeros(4,1);
    figure; hold on
    mismatch = false(1,size(split_vals,1));
    for k = 1:3
        h(k) = scatter(split_vals(actual_type==k,1), split_vals(actual_type==k,4)*16.5*16.5*23/10^9, '*', 'MarkerEdgeColor', C.colormap(2+k,:));
        mismatch = mismatch | (actual_type==k & pred_type ~=k);
        
    end
    
    h(4) = scatter(split_vals(mismatch,1), split_vals(mismatch,4)*16.5*16.5*23/10^9, 'o', 'MarkerEdgeColor', [0 0 0], 'sizeData', 72);
    
    prep_figure(gcf,gca, 'xlabel', '10th percentile depth (% IPL)', 'ylabel', 'Axonal arbor volume (um^3)', 'legend', {'BC3a', 'BC3b', 'BC4', 'Swapped'});
    set(gca, 'XTick', 10:5:35);
    
end
    
    