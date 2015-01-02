function [near_strat, far_strat] = stratification_by_depth_and_distance(target_cell, Q)
    
%     res = [16.5 16.5 25];
    
%     f = inline(' (30-108)/(2.7-1.5)/10^4.*(x-2.7*10^4)+30');
%     f = inline('-.006.*(x-2.7*10^4)+30');
    f = inline('1/(-2 - 1.2)/10^2*(x+2*10^4)+100');
%     f = inline('x');
%     f = inline('1/(5-11)/10^3*(x-11*10^5)');
%      f = inline('1/(1.1-3)/10^2*(x-3*10^4)');
     
%     f = inline('1/(5-10.5)/10^3*(x-10.5*10^5)');
        
    fn = ['./surface_points_trans/cell_' num2str(target_cell) '_surface.mat'];
%         if exist(fn,'file')
    load(fn);
    for k = 1:3
%         surface_points(:,k) = surface_points(:,k)*res(k);
    end                

    r = rand(size(surface_points,1),1)>.99;
    figure;        
    scatter(surface_points(r,2),surface_points(r,3),1);
    title(['jcell - ' num2str(target_cell)]);
    
    
%     surface_points = surface_points*Q';
    
    [min_z min_ind] = min(surface_points(:,1));
    prox_point = mean(surface_points(surface_points(:,1)==min_z,:),1);
    
    
    dist_from_prox = (surface_points(:,1) - prox_point(1)).^2 + ...
        (surface_points(:,2) - prox_point(2)).^2 + ...
        (surface_points(:,3) - prox_point(3)).^2;

    [max_dist max_ind] = max(dist_from_prox);
    dist_point = surface_points(max_ind,:);
    dist_vec = dist_point - prox_point;
    dist_vec = dist_vec(2:3) / sqrt(sum(dist_vec(2:3).^2));
    
    dist = (surface_points(:,2) - prox_point(2)).*dist_vec(1) + ...
        (surface_points(:,3) - prox_point(3)).*dist_vec(2);
    
    depth = f(surface_points(:,1));
%     depth_ticks = min(depth):max(depth);
%     depth = floor(depth-min(depth)+1);
%     
%     my_plot = zeros(max(depth),100);
%         
%     for n = 1:length(depth)
%         my_plot(depth(n),dist(n)) = my_plot(depth(n),dist(n))+1;
%     end
%     
%     figure;imagesc(my_plot); colorbar;
%     title(['number of surface points by depth and distance for ' num2str(10010)])
    
    
%     dist(depth<-19) = [];
%     depth(depth<-19) = [];
        
    x = -9:80;
    
    dist(depth>max(x)) = [];
    depth(depth>max(x)) = [];

    dist_split = median(dist);
    figure; subplot(2,1,1);
    near_strat = hist(depth(dist <= dist_split), x);
    hist(depth(dist <= dist_split), x);
    title(['stratification for distance from ' num2str(target_cell) ' cell body < ' num2str(dist_split)])
%     xlim([-20, 120]);
    
    
    subplot(2,1,2);
    far_strat = hist(depth(dist > dist_split), x);
    hist(depth(dist > dist_split), x);
    title(['stratification for distance from ' num2str(target_cell) ' cell body > ' num2str(dist_split)])
%     xlim([-20, 120])
    
    
%         dist_ticks = min(dist):(max(dist)-min(dist))/99:max(dist);

    dist = dist - min(dist);
    dist = floor(99*dist/max(dist)+1);
    

    depth_ticks = min(depth):max(depth);
    depth = floor(depth-min(depth)+1);
    
    my_plot = zeros(max(depth),100);
        
    for n = 1:length(depth)
        my_plot(depth(n),dist(n)) = my_plot(depth(n),dist(n))+1;
    end
    
    figure;imagesc(my_plot); colorbar;
    title(['number of surface points by depth and distance for ' num2str(10010)])
    
    
end
    
    
    