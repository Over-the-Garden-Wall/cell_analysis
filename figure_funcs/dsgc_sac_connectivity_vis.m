function dsgc_sac_connectivity_vis(dsgc, sac_nums, res)

    
    sac_sust_range = [0 40] * 1000 / res;
    sac_trans_range = [40 80] * 1000 / res;    
    conv_filter_size = 11;
    
    conts = cell(length(sac_nums),1);
    soma_locs = zeros(length(sac_nums), 3);
    
    
    for n = 1:length(sac_nums);
        c = detect_vericose_contacts(sac_nums(n), 500, 200, 50000);
        conts{n} = c(2:end, c(1,:)==dsgc);
        
        c_d = cell_data(sac_nums(n));
        soma_locs(n,:) = c_d.get_midpoint;
    end
    
    vol_min = min(soma_locs)-1;
    vol_max = max(soma_locs)+1;
    
    im_sz = ceil((vol_max(2:3)-vol_min(2:3))/res);
    
    exc = zeros(im_sz);
    sus_inh = zeros(im_sz);
    trans_inh = zeros(im_sz);
    
    for n = 1:length(sac_nums)
        soma_adj_pos = ceil((soma_locs(n,:) - vol_min)/res);
        
        for k = 1:size(conts{n},2)
            cont_adj_pos = ceil((conts{n}(2:4,k)' - vol_min)/res);
            theta = atan2(cont_adj_pos(3)-soma_adj_pos(3), cont_adj_pos(2)-soma_adj_pos(2));
            d = round(sqrt((cont_adj_pos(3)-soma_adj_pos(3))^2 + (cont_adj_pos(2)-soma_adj_pos(2))^2));
            
            for t = 0:d
                loc = round([soma_adj_pos(2) + t*cos(theta), soma_adj_pos(3) + t*sin(theta)]);
                if t >= sac_sust_range(1) && t <= sac_sust_range(2)
                    sus_inh(loc(1), loc(2)) = sus_inh(loc(1), loc(2)) + 1;
                end
                if t >= sac_trans_range(1) && t <= sac_trans_range(2)
                    trans_inh(loc(1), loc(2)) = sus_inh(loc(1), loc(2)) + 1;
                end
            end
        end
    end
       
    c_d = cell_data(dsgc);
    p = c_d.get_surface;
    p(:,2) = ceil((p(:,2) - vol_min(2))/res);
    p(:,3) = ceil((p(:,3) - vol_min(2))/res);
    
    for n = 1:size(p,1);
        exc(p(n,2), p(n,3)) = 1;
    end
    
    
    %make convolutional_filer
    g = gausswin(conv_filter_size);
    K = g*g';
    K = K/sum(K(:));
    
    ims = cell(4,1);
    
    ims{3} = conv2(exc, K);
    ims{2} = conv2(sus_inh, K);
    ims{1} = conv2(trans_inh, K);
    
    for n = 1:3
        ims{n} = ims{n} / max(ims{n}(:));
    end
    
    ims{4} = cat(3, ims{1}, ims{2}, ims{3});
    
    
    im_com = zeros(3,2);    
    for n = 1:3
        x_sum = 0;
        y_sum = 0;
        for y = 1:size(ims{n},1)
            for x = 1:size(ims{n},2);
                x_sum = x_sum + x*ims{n}(y,x);
                y_sum = y_sum + y*ims{n}(y,x);
            end
        end
        x_sum = x_sum / sum(ims{n}(:));
        y_sum = y_sum / sum(ims{n}(:));
        
        im_com(n,:) = [y_sum, x_sum];
    end    
    
    close all
    
    figure;
    titles = {'transient inhibitory RF', 'sustained inhibitory RF', 'excitatory RF', 'RGB combined'};
    
    for n = 1:4
        subplot(2, 2, n); 
        imagesc(ims{n}); hold all
        if n < 4
            scatter(im_com(n,2), im_com(n,1), 50, '*', 'markerEdgeColor', [0 0 0]);
        else
            scatter(im_com(1,2), im_com(1,1), 50, '*', 'markerEdgeColor', [0 1 1]);
            scatter(im_com(2,2), im_com(2,1), 50, '*', 'markerEdgeColor', [1 0 1]);
            scatter(im_com(3,2), im_com(3,1), 50, '*', 'markerEdgeColor', [1 1 0]);            
        end
        title(titles{n}); set(gca, 'XDir', 'reverse')
    end
    
    
    
    
end
    
            
            
    
    
        