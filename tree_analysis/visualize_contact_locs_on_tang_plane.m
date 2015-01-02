function visualize_contact_locs_on_tang_plane(cell_nums, cont_types, kernel_size)
    C = get_constants;

    MAX_SIZE = 1000;
    
    %assume max size of 1000 or so
    num_types = length(cont_types);
    cont_map = zeros(2*MAX_SIZE+1,2*MAX_SIZE+1, num_types);
    
    for n = cell_nums
        cell_dat = cell_data(n);
        
        cell_mid = cell_dat.get_midpoint(true);
        
        c = double(cell_dat.contacts);
        type_num = zeros(1,size(c,2));
        for k = 1:num_types
            for t = 1:size(c,2)
                if any(c(t)==C.type.(cont_types{k}))
                    type_num(t) = k;
                end
            end
        end
        
        for d = 1:2
            c(d+3,:) = round((c(d+3,:) - cell_mid(d+1))/1000);
        end
        
        for k = 1:num_types
            for t = find(type_num==k)
                cont_map(MAX_SIZE+c(4,t)+1,MAX_SIZE+c(5,t)+1,k) = ...
                    cont_map(MAX_SIZE+c(4,t)+1,MAX_SIZE+c(5,t)+1,k) + c(2,t);
            end
        end
    end
    
    K = gausswin(kernel_size);
    K = K*K';
    
    for k = 1:num_types
        cont_map(:,:,k) = conv2(cont_map(:,:,k), K, 'same');
    end
    
    total_map = sum(cont_map,3);
    lims = zeros(2);
    for d = 1:2
        t = sum(total_map,d);
        lims(d,1) = find(t,1,'first');
        lims(d,2) = find(t,1,'last');
    end
    lim = max(abs(lims(:)-MAX_SIZE-1));
    
    cont_map = cont_map((-lim:lim)+MAX_SIZE+1, (-lim:lim)+MAX_SIZE+1,:);
    
    c = colormap('Lines');
    
    for k = 1:num_types;
        figure;
        imagesc(cont_map(:,:,k), [0 max(cont_map(:))]);
        colorbar;
        title(cont_types{k});
        
        if length(cell_nums) == 1
            p = cell_dat.get_surface;
            p = p(:,2:3);
            for d = 1:2
                p(:,d) = p(:,d) - cell_mid(d+1);
            end
%             p(:,2) = -p(:,2);
            p = round(p/1000 + lim + 1);           
            p = unique(p,'rows');
            
            hold on;
            scatter(p(:,2),p(:,1),1,'markerEdgeColor', [1 1 1]);
        end
           
    end
    
    
    colored_im = zeros(size(cont_map,1), size(cont_map,2), 3);
    for d = 1:3
        for k = 1:num_types
            colored_im(:,:,d) = colored_im(:,:,d) + c(k,d)*cont_map(:,:,k);
        end
    end
    
    figure; imagesc(colored_im/max(colored_im(:)));
end