function [x y y2] = get_distribution_by_distance(target_cell, contacting_cells, R, Q)

res = [16.5 16.5 25];



    fn = ['./surface_points/cell_' num2str(target_cell) '_surface.mat'];
%         if exist(fn,'file')
    load(fn);
    for k = 1:3
        surface_points(:,k) = surface_points(:,k)*res(k);
    end                

    surface_points = surface_points*Q';
    
    [min_z min_ind] = min(surface_points(:,1));
    prox_point = mean(surface_points(surface_points(:,1)==min_z,:),1);
    dist_from_prox = (surface_points(:,1) - prox_point(1)).^2 + ...
        (surface_points(:,2) - prox_point(2)).^2 + ...
        (surface_points(:,3) - prox_point(3)).^2;

    [max_dist max_ind] = max(dist_from_prox);
    dist_point = surface_points(max_ind,:);
    dist_vec = dist_point - prox_point;
    
    
    




new_R = double(R);
for k=1:3
    new_R(3+k,:) = new_R(3+k,:)*res(k);
end
new_R(4:6,:) = Q*new_R(4:6,:);


is_good = false(1,size(R,2));    
for k = contacting_cells; 
        is_good = is_good | (R(1,:)==target_cell & R(2,:)==k);
        is_good = is_good | (R(2,:)==target_cell & R(1,:)==k);
        
end

new_R = new_R(:,is_good)';
R_dist = (new_R(:,5) - prox_point(2))*dist_vec(2) + (new_R(:,6) - prox_point(3))*dist_vec(3);
R_dist = R_dist/max(R_dist)*100;
R_area = new_R(:,3);

[R_dist sort_inds] = sort(R_dist);
R_area = R_area(sort_inds);

x = floor(R_dist(1)):ceil(R_dist(end));
R_dist = round(R_dist);


%consolidate R_area
for n = length(R_dist)-1:-1:1
    if R_dist(n) == R_dist(n+1)
        R_dist(n+1) = [];
        R_area(n) = R_area(n)+R_area(n+1);
        R_area(n+1) = [];
    end
end


y = zeros(1,length(x));

y(R_dist-x(1)+1) = R_area;

y2 = conv(y, gausswin(20)/sum(gausswin(20)), 'same');

% x = x -;

% f = inline('-.006*(x-2.7*10^4)+30');
% x = f(x);


figure;subplot(1,2,1); plot(x,y2)
title('smoothed contact area by distance')
xlabel('distance from soma on proximal-distal axis')
subplot(1,2,2); plot(x,y)
title('non-smoothed contact area by distance')


end