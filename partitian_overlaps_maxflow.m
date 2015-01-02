function colors = partitian_overlaps_maxflow(M)
    
    M = M-min(M(:)); %must be non-negative
    [x, y] = meshgrid(1:size(M,1), 1:size(M,2));
    


    colors = zeros(size(M,1),1);
    colors(1) = 1;
    
    flow_val = zeros(size(M,1),1);
    flow_val(1) = Inf;
    for k = 2:size(M,1);
        tM = M;
        tM([k end],:) = tM([end k],:);
        tM(:,[k end]) = tM(:, [end k]);
        
        tM(y>=x) = 0;
        tM = sparse(tM);
        
        flow_val(k) = graphmaxflow(tM, 1, size(M,1));
    end
    [dummy, k] = min(flow_val);
    
    tM = M;
    tM([k end],:) = tM([end k],:);
    tM(:,[k end]) = tM(:, [end k]);
        
    tM(y>=x) = 0;
    tM = sparse(tM);
        
    [dummy, dummy2] = graphmaxflow(tM, 1, size(M,1));
end