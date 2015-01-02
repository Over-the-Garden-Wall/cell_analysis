function visualize_overlaps(cell_nos, rad)
    num_cells = length(cell_nos);
    hull_points = cell(num_cells,1);

    intersection_levels = cell(1,1);
    for k = 1:length(cell_nos)
        c_d = cell_data(cell_nos(k));
%         p = c_d.get_surface;
%         p = p(1:100:end,2:3);
%         hull_points{k} = make_locally_convex_hull(p, rad, true);
        [hull_points{k}(:,1), hull_points{k}(:,2)] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
%         h = [];
%         [h(:,1), h(:,2)] = polybool('union', hull_points{k}(:,1), hull_points{k}(:,2), intersection_levels{1}(:,1), 
    end
    
    
    intersection_levels{1} = hull_points;
    ilev = 2; 
%     while 1
    for ilev = 2:4
        k = 1;
        for m = 1:length(intersection_levels{ilev-1})
            h1 = intersection_levels{ilev-1}{m};
            for n = m+1:length(intersection_levels{ilev-1})
                h2 = intersection_levels{ilev-1}{n};
                h = polybool('intersection', h1(:,1), h1(:,2), h2(:,1), h2(:,2));
                if ~isempty(h);
                    h = [];
                    [h(:,1), h(:,2)] = polybool('intersection', h1(:,1), h1(:,2), h2(:,1), h2(:,2));
                    intersection_levels{ilev}{k} = h;
                    k = k+1;
                end
            end
        end
        if k == 1
            break
        end
        
        is_not_redundant = true(length(intersection_levels{ilev}),1);
        for m = 1:length(intersection_levels{ilev})
            um = unique(intersection_levels{ilev}{m}(:));
            for n = m+1:length(intersection_levels{ilev})
                un = unique(intersection_levels{ilev}{n}(:));
                if length(um)==length(un) && all(um==un)
                    is_not_redundant(n) = false;
                end
            end
        end
        intersection_levels{ilev} = intersection_levels{ilev}(is_not_redundant);
        
%         ilev = ilev + 1;
    end
        
    c = colormap('jet');
    figure; hold all
    for k = 1:ilev-1;
        for m = 1:length(intersection_levels{k})
            fill(intersection_levels{k}{m}(:,1), intersection_levels{k}{m}(:,2), c(round(k/ilev*size(c,1)),:), 'lineStyle', 'none');
        end
    end
    
    for k = 1:length(cell_nos)
        text( mean(hull_points{k}(:,1)), mean(hull_points{k}(:,2)), num2str(cell_nos(k)));
    end
        
end