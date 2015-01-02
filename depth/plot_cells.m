function plot_cells(cell_nums, proj_plane, density, c, use_ipl_depth)

    if ~exist('use_ipl_depth','var') || isempty(use_ipl_depth)
        use_ipl_depth = false;
    end

%     res = [16.5 16.5 25];
    C = get_constants;


    num_cells = length(cell_nums);

    proj_axes = 1:3;
    proj_axes(proj_plane) = [];
    
%     figure;
    hold on
%     c = colormap('Lines');
    
   
    for n = 1:num_cells
        
        fn = [C.point_dir '/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            
            if use_ipl_depth
                surface_points(:,1) = C.f(surface_points(:,1));
            end
%             surface_points(:,1) = C.f(surface_points(:,1));
%             disp(['loaded ' num2str(size(surface_points,1)) ' points from cell ' num2str(cell_nums(n))])
            
            
            r = rand(size(surface_points,1),1) < density;
            
            scatter(surface_points(r,proj_axes(1)), surface_points(r,proj_axes(2)), 3, c);
            
        end
    end
end