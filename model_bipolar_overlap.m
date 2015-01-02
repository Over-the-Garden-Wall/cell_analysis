function ol_mat = model_bipolar_overlap
    
    C = get_constants;
    
    Jn = C.type.j;
    
    ct = {'t1','t2','t3a','t3b','t4'};
    
    %density per mm^2
    t_density = [2233 3212 1866 3254 3005];
%     t_sd = [399 543 605 536 573];
    
    num_types = length(ct);
    
    %first, make a density graph of the J cells. bin by um^3
    j_dens = zeros(100, 300,300); %z in IPL depth
    for n = Jn
        cell_dat = cell_data(n);
        p = cell_dat.get_surface;
        
        p(:,1) = round(C.f(p(:,1)));
        is_valid = p(:,1) >= 1 & p(:,1) <= 100;
        p = p(is_valid,:);
        
        mid_point = cell_dat.get_midpoint(true);
        
        p(:,2) = round((p(:,2) - mid_point(2))/1000)+151;
        p(:,3) = round((p(:,3) - mid_point(3))/1000)+151;
        
        
        for k = 1:size(p,1)
            j_dens(p(k,1),p(k,2),p(k,3)) = j_dens(p(k,1),p(k,2),p(k,3))+1;
        end       
    end
    j_dens = j_dens/length(Jn);
    
    t_dens = zeros(100,num_types);
    
    for type = 1:num_types;
        Tn = C.type.(ct{type});
        for n = Tn
            cell_dat = cell_data(n);
            p = cell_dat.get_surface;
        
            p(:,1) = round(C.f(p(:,1)));
            is_valid = p(:,1) >= 1 & p(:,1) <= 100;
            p = p(is_valid,:);
            
            t_dens(:,type) = t_dens(:,type) + hist(p(:,1),1:100)'/length(Tn);
        end
        t_dens(:,type) = t_dens(:,type)*t_density(type);
    end
    
    ol_mat = zeros(300, 300, num_types);
    
    for type = 1:num_types
        for k = 1:size(j_dens,1);
            ol_mat(:,:,type) = ol_mat(:,:,type) + squeeze(j_dens(k,:,:))*t_dens(k,type);
        end
    end
end
    
    
            