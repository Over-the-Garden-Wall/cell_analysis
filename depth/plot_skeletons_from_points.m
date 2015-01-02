function plot_cells(cell_nums, Q, proj_plane)

    res = [16.5 16.5 25];

    num_cells = length(cell_nums);

%     proj_axes = 1:3;
%     proj_axes(proj_plane) = [];
    
    figure; hold on
    c = colormap('Lines');
    

    f = inline(' -.006.*(x-2.7*10^4)+30');
    
    for n = 1:num_cells
        
        
        im = points2full(cell_nums(n), [3 3 3]);
        
        
        if ~isempty(im)
            
            [nodes edges] = skeletonize_by_grid(im, 32*ones(1,3));
            
            for k = 1:3
                nodes(:,k) = nodes(:,k) * res(k) * 3;
            end
            
            
            nodes = nodes*Q';
            
            nodes(:,1) = f(nodes(:,1));
%             disp(['loaded ' num2str(size(surface_points,1)) ' points from cell ' num2str(cell_nums(n))])
            
            plot_skeleton(nodes,edges,proj_plane, c(mod(n,size(c,1))+1,:));
            
            
        end
    end
end