function nearby_inds = find_nearby_points(query_p, p, min_dist)

    min_p = min(p);
    max_p = max(p);
    
    p_space = max_p-min_p;
    
    num_p = size(p,1);
    num_q = size(query_p,1);    
    num_d = size(p,2);
    
    
    best_num_bins = round(sqrt(num_q));
    
    
    space_fact = (best_num_bins/prod(p_space))^(1/num_d);
    
    num_bins = round(p_space*space_fact);
    num_bins = max([num_bins; ones(1,num_d)]);
    
    total_num_bins = prod(num_bins);
    bin_min = zeros(total_num_bins,3);
    bin_max = zeros(total_num_bins,3);
    
    eval_str = '[';
    for d = 1:num_d
        eval_str = [eval_str 'bin_list(:,' num2str(d) '), '];
    end
    
    eval_str = [eval_str(1:end-2) '] = ind2sub(num_bins, 1:total_num_bins);'];
    eval(eval_str);
    bin_spacing = p_space./num_bins;
    
    
    
    for n = 1:total_num_bins        
        bin_min(n,:) = (bin_list(n,:)-1).*bin_spacing + min_p;
        bin_max(n,:) = bin_list(n,:).*bin_spacing + min_p;
    end
    
    bin_min(bin_list==1) = -Inf;
    for d = 1:num_d
        bin_max(bin_list(:,d)==num_bins(d),d) = Inf;
    end
        
    sub_p = cell(total_num_bins,1);
    sub_p_inds = cell(total_num_bins,1);
    
    for n = 1:total_num_bins % O(num_p*total_num_bins)
        is_in = true(num_p,1);
        for d = 1:num_d
            is_in = is_in & p(:,d)>bin_min(n,d)-min_dist & p(:,d)<bin_max(n,d)+min_dist;
        end
        
        sub_p{n} = p(is_in,:);
        sub_p_inds{n} = find(is_in);
    end
    
    nearby_inds = cell(num_q,1);
    
    for n = 1:num_q
        is_in = true(total_num_bins,1);
        for d = 1:num_d
            is_in = is_in & query_p(n,d)>bin_min(:,d) & query_p(n,d)<bin_max(:,d);
        end
        my_bin = find(is_in,1,'first');
        
        dist = zeros(size(sub_p{my_bin},1),1);
        for d = 1:num_d
            dist = dist + (sub_p{my_bin}(:,d)-query_p(n,d)).^2;
        end
        dist = sqrt(dist);
        nearby_inds{n} = sub_p_inds{my_bin}(dist<min_dist);
    end

end
        
        
    