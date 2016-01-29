function figure_sac_sideview

    C = get_constants;
    
    hist_x = -10:100;    
    
    sac_type = {'sure_off_sac', 'on_sac'};
    plot_c = [0 1 1; 1 1 0];
    
    rel_d = [1 2];
    v_bounds = C.volume_bounds(:,rel_d);
    
    pixel_size = 100*ones(1,2);
    
    im_size = ceil((v_bounds(2,:) - v_bounds(1,:))./pixel_size) + 30; %30 wiggle room points
    im = zeros([im_size 2]);
    
    sampling_rate = [10 100];   
    
    hist_data = zeros(length(hist_x),2,2);
    
    load(C.trans_loc);    
    
    for align_c = 1:2
    
    %     figure; hold all;
        for l = 1:2

            sac_nums = C.type.(sac_type{l});

            for c = sac_nums
                c_d = cell_data(c);
                
                if align_c == 1
                    p = c_d.get_surface;
                    p = p(round(sampling_rate(align_c)/2:sampling_rate(align_c):end), rel_d);
                else
%                     p = c_d.get_volume(ones(3,1) * [-Inf Inf]);
                    load([C.surface_point_dir 'cell_' num2str(c) '_surface.mat']);
                    p = surface_points(round(sampling_rate(align_c)/2:sampling_rate(align_c):end), :);
                    for n = 1:size(p,2)
                        p(:,n) = p(:,n)*C.res(n);
                    end
                    p = p*T.global_Q';
                    p = p(:,rel_d);
                    p(:,1) = p(:,1) - mean(T.off_chat);
                end
                
                
                

                dpth = C.f(p(:,1));
                d_hist = hist(dpth, [hist_x(1)-1, hist_x, hist_x(end)+1]);
                hist_data(:,align_c,l) = hist_data(:,align_c,l)+d_hist(2:end-1)';
                
                
                for d = 1:2
                    p(:,d) = ceil((p(:,d) - v_bounds(1,d)) / pixel_size(d));
                end


                for n = 1:size(p,1)
                    im(p(n,1), p(n,2),l) = im(p(n,1), p(n,2),l) + 1;
                end

    %             s = load([C.skele_dir 's' num2str(c) '.mat']);
    %             plot_skeleton(s, 3, plot_c(l,:));
            end
        end

        im = log(im+1);

        im_max = [max(max(im(:,:,1))) max(max(im(:,:,2)))];
        disp_im = zeros([im_size, 3]);
        for d = 1:2
            for c = 1:3
                disp_im(:,:,c) = disp_im(:,:,c) + im(:,:,d)/im_max(d) * plot_c(d,c);
            end
        end
        disp_im(disp_im>1) = 1;
        figure; imagesc(1-disp_im)
        imwrite(1-disp_im, [C.image_dir '/sac_sideview_' num2str(align_c) '.png'], 'png', 'Transparency', [1 1 1]);

    end
    
    
%     pixel_size = 5*ones(1,2);
    for k = 1:4
        hist_data(:,k) = hist_data(:,k) / sum(hist_data(:,k));
    end
    figure; hold all
    plot(1:100, hist_data(11:110,[1 3]), 'LineWidth', 2);
    plot(1:100, hist_data(7:106,[2 4]), 'LineWidth', 2);
    prep_figure(gcf,gca, 'legend', {'Aligned OFF', 'Aligned ON', 'Raw OFF', 'Raw ON'});
    
    whitebg('black')
    set(gcf, 'Color', [0 0 0])
    set(gcf, 'InvertHardCopy', 'off');
    
    
    
end