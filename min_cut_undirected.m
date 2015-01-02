function cut = min_cut_undirected(M, src)
    
    
    num_nodes = size(M,1);

    cut_hist = cell(num_nodes,1);
    cut_vals = zeros(num_nodes,1);
    
    node_list = 1:num_nodes;
    
    cut_hist{1} = src;
    node_list(node_list==src) = [];    
    M = M([src, 1:src-1, src+1:end], [src, 1:src-1, src+1:end]);
    cut_vals(1) = sum(M(1,:));    
    M(1,1) = 0;
    
    
    
    for k = 2:num_nodes
        G = M;
        [x, y] = meshgrid(1:size(M,1), 1:size(M,2));
        G(y>=x) = 0;
        G = sparse(G);
        [dummy, worst_node] = min(M(1,2:end));
        flow_val = graphmaxflow(G, 1, worst_node+1);
        
        [dummy, tightest_node] = max(M(1,:));
        cut_hist{k} = [cut_hist{k-1} node_list(tightest_node-1)];
        M(:,1) = M(:,1)+M(:,tightest_node);
        M(1,:) = M(1,:)+M(tightest_node,:);
        M(1,1) = 0;
        M(tightest_node,:) = [];
        M(:,tightest_node) = [];
        node_list(tightest_node-1) = [];
        
        cut_vals(k) = sum(M(1,:));
        disp([cut_vals(k) flow_val])
    end
    
    
end
        
        

    

