function [s num_points] = strat_notrunk(cell_num)

    C = get_constants;    

    c_d = cell_data(cell_num);
    
    s = c_d.stratification * c_d.SA;
    old_s = s;
    
%     smooth_s = conv(s, ones(1,11)/11, 'same');
%     [a, def_not_trunk] = max(smooth_s);
%         
%     always_allowed = smooth_s > 2500;
%     
%     trunk_end = find(smooth_s(1:def_not_trunk) < 1000, 1, 'last');
%     
    
%     s((1:length(s))' <= trunk_end & ~always_allowed) = 0;

    skele = load([C.skele_dir 's' num2str(cell_num) '.mat']);
    sort_edges = sort(skele.edges(:));
    is_branch = false(size(skele.nodes,1), 1);
    for n = 1:length(sort_edges)-2
        if sort_edges(n) == sort_edges(n+2)
            is_branch(sort_edges(n)) = true;
        end
    end
    node_d = C.f(skele.nodes(is_branch,1));
    node_d(node_d < 5) = [];
    min_d = round(min(node_d)-5);    

    s(1:min_d - C.strat_x(1) + 1) = 0;
    figure; plot([old_s s]);
    num_points = sum(s);
    s = s/num_points;
    
end
