function dsgc_rf(gc)
    
    C = get_constants;
    c_d = cell_data(gc);
    
    h = cell(1,2);
    [h{:}] = poly2cw(c_d.hull_2d(:,1), c_d.hull_2d(:,2));
    
    ex_h = expand_hull(h, 10000);
    
    figure; hold all
    s = load([C.skele_dir 's' num2str(gc) '.mat']);
    
    plot([s.nodes(s.edges(:,1),2) s.nodes(s.edges(:,2),2)]', ...
        [s.nodes(s.edges(:,1),3) s.nodes(s.edges(:,2),3)]', ...
        'Color', .5*ones(1,3), 'lineWidth', 2);
    
    plot(ex_h{1}, ex_h{2}, 'k', 'lineWidth', 2);
    
    t3aex_h = spread_hull(h, 120:240, 50000, 30000);
    t2ex_h = spread_hull(h, 120:240, 90000, 30000);
    
    plot(t3aex_h{1}, t3aex_h{2}, 'r', 'lineWidth', 2);
    plot(t2ex_h{1}, t2ex_h{2}, 'b', 'lineWidth', 2);
end

function h = expand_hull(h, dist)

    p = zeros([360 length(h{1}) 2]);
    
    for n = 1:length(h{1})    
        p(:, n, 1) = cos((1:360)/180*pi)*dist + h{1}(n);
        p(:, n, 2) = sin((1:360)/180*pi)*dist + h{2}(n);
    end
    
    p = reshape(p, [size(p,1)*size(p,2), 2]);
    is_hull = convhull(p(:,1), p(:,2));
    
    [h{:}] = poly2cw(p(is_hull,1), p(is_hull,2));
end


function h = spread_hull(h, angles, dist, stretch)

    p = zeros([length(angles), 2, length(h{1}) 2]);
    
    for n = 1:length(h{1})    
        p(:, 1, n, 1) = cos(angles/180*pi)*(dist-stretch(1)/2) + h{1}(n);
        p(:, 1, n, 2) = sin(angles/180*pi)*(dist-stretch(1)/2) + h{2}(n);
        p(:, 2, n, 1) = cos(angles/180*pi)*(dist+stretch(1)/2) + h{1}(n);
        p(:, 2, n, 2) = sin(angles/180*pi)*(dist+stretch(1)/2) + h{2}(n);
    end
    
    p = reshape(p, [size(p,1)*size(p,2)*size(p,3), 2]);
    is_hull = convhull(p(:,1), p(:,2));
%     is_hull = true(size(p,1),1);
    
    [h{:}] = poly2cw(p(is_hull,1), p(is_hull,2));
end