function [scatter_dat] = connectivity_scatter(ref_cells, comp_cells)

    nc1 = length(ref_cells);
    total_contact = zeros(nc1,length(comp_cells));
    total_overlap = zeros(nc1,length(comp_cells));
    
    for t = 1:length(comp_cells)
        nc2 = length(comp_cells{t});

        overlap_matrix = zeros(nc1, nc2);
        connectivity_matrix = zeros(nc1, nc2);


        for m = 1:nc1
            c_d1 = cell_data(ref_cells(m));
            my_conns = double(c_d1.contacts);
            for n = 1:nc2
                c_d2 = cell_data(comp_cells{t}(n));

                h1 = c_d1.hull_2d;
                [h1(:,1) h1(:,2)] = poly2cw(h1(:,1), h1(:,2));

                h2 = c_d2.hull_2d;
                [h2(:,1) h2(:,2)] = poly2cw(h2(:,1), h2(:,2));

                h = [];
                try
                    [h(:,1), h(:,2)] = polybool('intersection', h1(:,1), h1(:,2), h2(:,1), h2(:,2));
                    ol = polyarea(h(:,1), h(:,2));
                catch ME
                    ol = 0;
                end

                conn_is_c2 = my_conns(1,:) == comp_cells{t}(n);

                total_overlap(m,t) = total_overlap(m,t)+ol;
                total_contact(m,t) = total_contact(m,t)+sum(my_conns(2,conn_is_c2));                
            end
        end

    end
    
    scatter_dat = total_contact./total_overlap;
%     for t = 2:size(scatter_dat,2)
%         figure; scatter(scatter_dat(:,1), scatter_dat(:,t));
%     end
    figure; hold all
%     for t = 1:size(scatter_dat,2)
%         scatter(t*ones(size(scatter_dat,1),1), scatter_dat(:,t), '*');
%     end
leg_lbl = cell(size(scatter_dat,1),1);
    
    for t = 1:size(scatter_dat,1);
        scatter(1:size(scatter_dat,2), scatter_dat(t,:), '*');
        leg_lbl{t} = num2str(ref_cells(t));
    end
    legend(leg_lbl);
    
    
end
    