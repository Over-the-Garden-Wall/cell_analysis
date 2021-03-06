function type_gallery(types, type_labels)
    
    P_SPARSITY = 500;
    CELL_ID_MAX = 100000;

    LINE_MAX = 10^5;
    
    close all
    
    C = get_constants;
    num_types = length(types);
    
    props_to_copy = {'XLim', 'YLim', 'XTick', 'YTick', 'XTickLabel', 'YTickLabel'};
    
    %prepare matrices
    
    type_nums = cell(1,num_types);
    for k = 1:num_types
        type_nums{k} = C.type.(types{k});
    end
        
    num2type = zeros(CELL_ID_MAX, 1);
    for k = 1:num_types
        num2type(C.type.(types{k})) = k;
    end
            
    M = zeros([num_types num_types 2]); %pix1 in overlap, pix2 in overlap
    M(:,:,1) = points_within_hull_intersect(type_nums, type_nums);
    M(:,:,2) = M(:,:,1)';
    
    total_contact = contact_total(type_nums, type_nums);
    
    %begin making panels
    
    for k = 1:num_types
        cns = C.type.(types{k});
        colmap = make_colormap(length(cns), 5);
        
%         skeleton panels
        for proj_plane = 1:3
            if proj_plane == 1
                figure; h = subplot(2,1,1); hold all; 
                h2 = subplot(2,1,2); hold all;
                axes(h)
                
            elseif proj_plane == 2
                axes(h2)
            else
                figure; hold all;
                prep_figure(gcf,gca);
            
            end
                warning('single skeleton hack in place');
            for c = 1:0 %length(cns)
                s = load([C.skele_dir 's' num2str(cns(c)) '.mat']);
                
                if size(s.edges, 1) * length(cns) > LINE_MAX
                    s = reduce_skeleton(s, LINE_MAX / (size(s.edges, 1) * length(cns)));
                end
                
                s.nodes = s.nodes(:, [2, 3, 1]);
                s.nodes(:,3) = C.f(s.nodes(:,3));
                s.nodes(:,1:2) = s.nodes(:,1:2)/1000;
%                 s.nodes(:,1) = 300 - s.nodes(:,1);
                plot_skeleton(s, proj_plane, colmap(c,:));                       
            end
            
            
            if proj_plane ~= 3
                set(gca, 'YLim', [-20 120], 'YTick', [0 100], 'YTickLabel', {'1', '0'});                
                fig_sz = [1280 640];
                ax_sz = [120 19]*10;
                rel_sz = ax_sz ./ fig_sz;
%                 rel_sz(2) = rel_sz(2)/2;
                
                vert_space = (1 - rel_sz(2)*2) / 3;
                hor_space = (1 - rel_sz(1))/2;                                
                
                set(gcf, 'position', [0 0 fig_sz]);
                if proj_plane == 2
                    set(gca, 'position', [hor_space, vert_space, rel_sz]);
                else
                    set(gca, 'position', [hor_space, 1-vert_space-rel_sz(2), rel_sz]);
                end
                plot([0 300], [0 0] , '--', 'Color', [0 0 0]);
                plot([0 300], [100 100], '--', 'Color', [0 0 0]);
                plot([0 300], [28 28], ':', 'Color', [0 0 0]);
                plot([0 300], [62 62], ':', 'Color', [0 0 0]);

                plot([10 60], [115 115], 'LineWidth', 3, 'Color', [0 0 0]);

                %[0 0 120 19]*10);
            else
%                 set(gca, 'YLim', [0 3]*10^2, 'YTick', (0:3)*10^2);
                set(gca, 'YLim', [0 3]*10^2, 'YTick', []);
                plot([10 60], [290 290], 'LineWidth', 3, 'Color', [0 0 0]);
                
                
            end
%             set(gca, 'XLim', [0 3]*10^2, 'XTick', (0:3)*10^2);
            set(gca, 'XLim', [0 3]*10^2, 'XTick', []);
            
            
            
            warning('rasterization off');
%             rasterize_figure(gca);
            
            if proj_plane == 3
                xl = get(gca, 'XLim');
                yl = get(gca, 'YLim');
                
                plot([1 xl(2)], [1 1] , 'Color', [0 0 0]);
                plot([1 xl(2)], [yl(2) yl(2)] , 'Color', [0 0 0]);
                plot([1 1], [1 yl(2)] , 'Color', [0 0 0]);
                plot([xl(2) xl(2)], [1 yl(2)], 'Color', [0 0 0]);
            end
            
                
            
            set(gca, 'fontSize', 20);
            
            if proj_plane > 1
                save_fig([C.image_dir 'type_gallery_' types{k} '_skeles_' num2str(proj_plane) '.eps']);
            end
%             saveas(gcf, [C.image_dir 'type_gallery_' types{k} '_skeles_' num2str(proj_plane) '.png'], 'png');
%             close all
        end 
        
        
        %mosaic panel
        visualize_overlaps(cns, []);        
        prep_figure(gcf,gca);
%         set(gca, 'YLim', [0 3]*10^5, 'YTick', (0:3)*10^5, 'YTickLabel', {'0', '100', '200', '300'});
%         set(gca, 'XLim', [0 3]*10^5, 'XTick', (0:3)*10^5, 'xTickLabel', {'0', '100', '200', '300'});
        set(gca, 'YLim', [0 3]*10^5, 'YTick', []);
        set(gca, 'XLim', [0 3]*10^5, 'XTick', []);

        plot([1 300000], [1 1]*500 , 'Color', [0 0 0]);
            plot([1 300000], [300000 300000]*.999 , 'Color', [0 0 0]);
            plot([1 1], [1 300000] , 'Color', [0 0 0]);
            plot([300000 300000], [1 300000], 'Color', [0 0 0]);
        
            set(gca, 'YDir', 'reverse');
            
        save_fig([C.image_dir 'type_gallery_' types{k} '_mosaic.eps']);
%         saveas(gcf, [C.image_dir 'type_gallery_' types{k} '_mosaic.png'], 'png');
        
        %connectivity panels
        figure;             
        prep_figure(gcf,gca);
            
        for n = 1:2
            subplot(2,1,n);
            bar(3.5 + log(total_contact(k,:)./M(k,:,n))/log(10),1);
            set(gca, 'YLim', [0 3], 'YTick', .5:2.5, 'YTickLabel', {'.1', '1', '10'});                            
            set(gca, 'XTick', 1:100, 'XLim', [.5 size(total_contact,2) + .5]);%1:num_types, 'XTickLabel', type_labels);
            set(gca, 'fontSize', 20);
%             set(gcf, 'position', [0 0 120 19]*10);
%             temp_axis = gca; 
%             plot_handles = get(temp_axis, 'Children');
%             copyobj(plot_handles, ax_handle(4+n));
%             for c = 1:length(props_to_copy)
%                 set(ax_handle(4+n), props_to_copy{c}, get(temp_axis, props_to_copy{c}));
%             end
            
        end        
        set(gcf, 'position', [0 0 fig_sz]);
        save_fig([C.image_dir 'type_gallery_' types{k} '_connectivity.eps']);
%             saveas(gcf, [C.image_dir 'type_gallery_' types{k} '_connectivity_' num2str(n) '.png'], 'png');        
%         end
    
        %stratification panel
%         plot_mean_strat({C.type.(types{k})}, [-10 110]); hold all;
        figure; hold all;
        h_range = -10:110;
        mean_plot = zeros(size(h_range));    
        
        for cn = 1:length(C.type.(types{k}))
            c = C.type.(types{k})(cn);
            c_d = cell_data(c);
            p = c_d.get_surface;
            d = C.f(p(:,1));
            d(d<h_range(1) | d>h_range(end)) = [];
            figure; plot_data = hist(d, h_range); close(gcf);            
            plot_data = plot_data/sum(plot_data);            
            plot(h_range, plot_data, 'Color', colmap(cn,:));
            mean_plot = mean_plot+plot_data;
        end
        mean_plot = mean_plot / sum(mean_plot);
        plot(h_range, mean_plot, 'Color', [0 0 0], 'lineWidth', 2);            
        
        prep_figure(gcf,gca);
        set(gca, 'YLim', [0 .2], 'YTick', [], 'YTickLabel', {' '});
        set(gca, 'XLim', [-10 110], 'XTick', 0:25:100);    
        save_fig([C.image_dir 'type_gallery_' types{k} '_strat.eps']);
%         saveas(gcf, [C.image_dir 'type_gallery_' types{k} '_strat.png'], 'png');                
%    
        figure; hold on;
        for cn = 1:length(C.type.(types{k}))
            c = C.type.(types{k})(cn);
            c_d = cell_data(c);
            p = c_d.get_surface;
            d = C.f(p(:,1));
            d(d<0) = [];
            d(d>100) = [];
            CA_p = unique(round(p(:,2:3)/1000), 'rows');
            CA = log(size(CA_p,1))/log(10);
            d = sort(d);
            r = [d(round(length(d)*.25)), d(round(length(d)*.75))];
            
            plot([r(1), r(2)], [CA CA], 'Color', colmap(cn,:));
%             h = scatter(d(round(length(d)*.25)), d(round(length(d)*.75)), 30);
%             set(h, 'markerEdgeColor', [0 0 0], 'markerFaceColor', [0 0 0]);
        end
        
        xtick = [0:1:4];
        xticklbl = cell(1,length(xtick));
        for n = 1:length(xticklbl)
            xticklbl{n} = ['10^' num2str(n-1)];
        end
        xticklbl{1} = '1';
        xticklbl{2} = '10';
        
        
        prep_figure(gcf, gca);
        set(gca, 'XLim', [0 100], 'XTick', [0:25:100]);                            
        set(gca, 'YLim', [0 xtick(end)], 'YTick', xtick, 'YTickLabel', xticklbl); 
        
        save_fig([C.image_dir 'type_gallery_' types{k} '_scatter.eps']);
%         
        set(gca, 'fontSize', 20);
            close all
% %     while gcf ~= main_fig
% %         close(gcf);
% %     end
% %     save_fig([C.image_dir 'type_gallery_' types{k} '_fullpage.eps'])
%         
% %         save_fig([C.image_dir 'type_gallery_' types{k} '_scatter.png']);
% %         saveas(gcf, [C.image_dir 'type_gallery_' types{k} '_scatter.png'], 'png');  
    end
    
end


        
        
            
            