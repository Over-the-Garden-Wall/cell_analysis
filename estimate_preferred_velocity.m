function activity = estimate_preferred_velocity(dist_domain, cell_prevalence, time_domain, type_activity)

    time_domain = time_domain*1000; %seconds -> milliseconds
    dist_domain = dist_domain*100; %microns -> hundredths of microns

    max_dist = max(dist_domain);

    speeds = 1:200; %hundredths of microns per millisecond
    
    start_point = ceil(dist_domain(1));
    time_of_sim = floor((max_dist-start_point)/speeds(1));
    time_of_activity = floor(max(time_domain)/speeds(1));
    
%     t_range = 1:time_of_sim;
    
    
    
    num_types = size(cell_prevalence,2);
    
    new_activity = zeros(time_of_activity,num_types);
    for tc = 1:time_of_activity
        new_activity(tc,:) = interp1q(time_domain,type_activity,tc);
    end
    
    new_prev = zeros(ceil(max(dist_domain)),num_types);
    for tc = start_point:size(new_prev,1);
        new_prev(tc,:) = interp1q(dist_domain,cell_prevalence,tc);
%         new_prev(tc,:) = new_prev(tc,:)/sum(new_prev(tc,:));
    
    end
    
    
%     for n = 1:num_types
% %         new_prev(:,n) = new_prev(:,n)/sum(new_prev(:,n));
%         new_activity(:,n) = new_activity(:,n)/sum(new_activity(:,n));
%     end        
    
    
    activity = zeros(time_of_sim+time_of_activity,length(speeds));
        
        
    for t = 1:time_of_sim
        time_range = t + (0:time_of_activity-1); 
        for s = 1:length(speeds);
            loc = ceil(t*speeds(s))+start_point;
            if loc <= max_dist
                for n = 1:num_types
                    activity(time_range,s) = activity(time_range,s) + new_prev(loc, n)*new_activity(:,n)*speeds(s);
                end
            end
        end
    end
    
                
            
end
    