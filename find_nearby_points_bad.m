function nearby_inds = find_nearby_points(query_p, p, min_dist)

    warning('this is slow, intersect is poorly written');

    num_d = size(p,2);
    px = cell(num_d,1);
    px_ind = cell(num_d,1);
    qp = cell(num_d,1);
    qp_ind = cell(num_d,1);
    min_ind = cell(num_d,1);
    max_ind = cell(num_d,1);
    
    num_p = size(p,1);
    num_q = size(query_p,1);
    
    nearby_inds = cell(num_q,1);
    
    for d = 1:num_d %this should be O(nlog(n)+mlog(m))
        
        [px{d} px_ind{d}] = sort(p(:,d));
        [qp{d} qp_ind{d}] = sort(query_p(:,d));
        [dummy qp_ind{d}] = sort(qp_ind{d});
        
        min_ind{d} = ones(length(qp{d}),1);
        max_ind{d} = zeros(length(qp{d}),1);                
        
    end
    
    for d = 1:num_d %this should be O(n+m)
        qk = 1;
        pk = 1;
        while 1
            if px{d}(pk) < qp{d}(qk)-min_dist
                pk = pk + 1;                
                if pk > num_p
                    min_ind{d}(qk:end) = pk;
                    break
                end
            else
                min_ind{d}(qk) = pk;                
                qk = qk + 1;
                if qk > num_q
                    break
                end
            end            
        end
        
        qk = 1;
        pk = 1;
        while 1
            if px{d}(pk) <= qp{d}(qk)+min_dist
                pk = pk + 1;                
                if pk > num_p
                    max_ind{d}(qk:end) = pk-1;
                    break
                end
            else
                max_ind{d}(qk) = pk-1;                
                qk = qk + 1;
                if qk > num_q
                    break
                end
            end            
        end
    end
            
    poss_inds = cell(num_d,1);
    for qk = 1:num_q % O(m * number_of_candidate_points)       
        for d = 1:num_d
            qd_ind = qp_ind{d}(qk);
            poss_inds{d} = px_ind{d}(min_ind{d}(qd_ind):max_ind{d}(qd_ind));
        end
        
        for d = 2:num_d
            poss_inds{1} = intersect(poss_inds{1},poss_inds{d});
        end
        
        cand_points = p(poss_inds{1},:);
        dist = zeros(size(cand_points,1),1);
        for d = 1:num_d
            dist = dist + (cand_points(:,d)-query_p(qk,d)).^2;
        end
        dist = sqrt(dist);
        
        nearby_inds{qk} = poss_inds{1}(dist <= min_dist);
    end
       
end
        
        
    