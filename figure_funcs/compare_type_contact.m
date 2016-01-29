function compare_type_contact(x_types, y_types, min_overlap)

    num_xtypes = length(x_types);
    num_ytypes = length(y_types);

    contDens_mean = zeros(num_xtypes, num_ytypes);
    contDens_ste = zeros(num_xtypes, num_ytypes);
    
    
    for xt = 1:num_xtypes        
        for yt = 1:num_ytypes
            [contact, overlap, pr_pred] = contactXoverlap(x_types{xt}, y_types{yt});            
            is_valid = overlap(:) > min_overlap;
            contact = contact(is_valid);
            overlap = overlap(is_valid);                        
            pr_pred = pr_pred(is_valid);
            
            contDens_mean(xt, yt) = mean((contact ./ overlap) - pr_pred);
            contDens_ste(xt, yt) = std((contact ./ overlap) - pr_pred) / sqrt(length(pr_pred)-1);
        end
    end
    figure;
    error_dot_plot(contDens_mean, contDens_ste);
    
end


function [cont, overlap, pr_pred] = contactXoverlap(xcells, ycells)
    C = get_constants;

    cont = zeros(length(xcells), length(ycells));
    overlap = zeros(length(xcells), length(ycells));
    pr_pred = zeros(length(xcells), length(ycells));

    for xc = 1:length(xcells)
        c_d = cell_data(xcells(xc));
        conts = double(c_d.contacts);
        h = cell(1,2);
        [h{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
        
        s = c_d.stratification;
        s = s(-C.strat_x(1):end);        
        
        Vdensity = c_d.V / polyarea(h{1}, h{2});
        
        for yc = 1:length(ycells)
            cont(xc,yc) = sum(conts(2, conts(1,:) == ycells(yc)));
            
            c_dy = cell_data(ycells(yc));
            hy = cell(1,2);
            [hy{:}] = poly2cw(c_dy.hull_2d(:,1), c_dy.hull_2d(:,2));
            
            h_intersect = cell(1,2);
            try
            [h_intersect{:}] = polybool('intersection', h{1}, h{2}, hy{1}, hy{2});
            
            overlap(xc, yc) = polyarea(h_intersect{1}, h_intersect{2});
            catch
            end
            
            Vdensityy = c_dy.V / polyarea(hy{1}, hy{2});
            
            sy = c_dy.stratification;
            sy = sy(-C.strat_x(1):end);
            
            strat_ol = sum(s(1:min(end, length(sy))).*sy(1:min(end, length(s))));
            pr_pred(xc, yc) = strat_ol * Vdensityy * Vdensity * 6;
            
        end
    end
end
                    