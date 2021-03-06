function path = skeletonize_sac(cell_num, dsmp_fact, separation_for_leaf, minimum_leaf_dist)
    tic
    
    C = get_constants;
    c = cell_data(cell_num);
    
    all_points = c.get_surface;
            
    p = unique(round(all_points/dsmp_fact),'rows');
    soma = c.get_midpoint(true);
    soma = soma/dsmp_fact;
    pmin = min(p);
    for d = 1:3;
        p(:,d) = p(:,d)-pmin(d)+1;
        soma(d) = soma(d) - pmin(d) + 1;
    end
    
    pmax = max(p);
    num_p = size(p,1);
    
    minimum_leaf_dist = minimum_leaf_dist/dsmp_fact;
    separation_for_leaf = separation_for_leaf/dsmp_fact;
    
    all_dist = euclid_dist(p, soma);    
    poss_leaf_p = p(all_dist>minimum_leaf_dist,:);
    [min_dist, S] = min(all_dist);
    
    p_graph = zeros(pmax+2);
    for n = 1:num_p
        p_graph(p(n,1)+1,p(n,2)+1,p(n,3)+1) = n;
    end
        
    [nhood(:,1) nhood(:,2) nhood(:,3)] = ind2sub([3 3 3],1:27);
    nhood = nhood - 1;
    nhood(all(nhood==1,2),:) = [];
    hoodSz = size(nhood,1);
    
    g_map = zeros(num_p,hoodSz);
    for n = 1:num_p
        for h = 1:hoodSz
            g_map(n,h) = p_graph(p(n,1)+nhood(h,1),p(n,2)+nhood(h,2),p(n,3)+nhood(h,3));
        end
    end
    
    [P hood_ind] = find(g_map);
    Q = g_map(g_map~=0);
    
    G = sparse(P,Q,ones(length(P),1));
    
    k = 1;
    while ~isempty(poss_leaf_p)
        dist = euclid_dist(poss_leaf_p, soma);
        max_dist = max(dist);
        
        T = find(all_dist==max_dist,1,'first');
        
        [d path_ind{k}] = graphshortestpath(G,S,T,'Directed', false, 'Method', 'BFS');       
        
        if isinf(d) %T is an orphan point
            dist = euclid_dist(poss_leaf_p, p(T,:));
            poss_leaf_p(dist <= separation_for_leaf,:) = [];
        else
            
            
        
            for n = 1:size(path_ind{k},2)
                dist = euclid_dist(poss_leaf_p, p(path_ind{k}(n),:));
                poss_leaf_p(dist <= separation_for_leaf,:) = [];
            end
        
        
%         path_ind_p = p(path_ind{k},:)+1;
        
%         figure; scatter(p(:,2)+1,p(:,3)+1,1)
%         
%         figure; scatter(p(:,2)+1,p(:,3)+1,1)
%         hold on; scatter(path_ind_p(:,2), path_ind_p(:,3), 1, 'r');
        
        
%             disp([num2str(k) ':   ' num2str(size(poss_leaf_p,1))])
            k = k+1;
        end
    end
    if isempty(path_ind{end})
        path_ind(end) = [];
    end
    
%     for k = 1:length(path_ind)
%         
%         
%     end
%         
%     warning('debug in place')
%     for n = 1:length(path_ind)
%         close all;
%         figure;
%         imagesc(squeeze(any(p_graph,1))); colormap gray;
%         hold on
%         scatter(p(path_ind{n},3), p(path_ind{n},2), '*r');
%         a=1;
%     end

    toc
    disp('converting to original coordinates and centering');
    tic
    
    for d = 1:3;
        p(:,d) = p(:,d)-1+pmin(d);
    end
    
    p = p*dsmp_fact;
      
    nearby_inds = find_nearby_points(p, all_points, 1000);
    
    for k = 1:size(p,1)
        p(k,:) = mean(all_points(nearby_inds{k},:));
    end
        
    path = cell(length(path_ind),1);
    for k = 1:length(path_ind)
        path{k} = p(path_ind{k},:);
    end
    
end
    

function dist = euclid_dist(p, ref)

    dist = sqrt((p(:,1)-ref(1)).^2 + ...
        (p(:,2)-ref(2)).^2 + ...
        (p(:,3)-ref(3)).^2);
    
end
    
    
    
    