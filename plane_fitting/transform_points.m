function transform_points(in_dir,out_dir,T, force_transform)

    if ~exist('force_transform','var') || isempty(force_transform)
        force_transform = false;
    end
    point_dirs = dir(in_dir);
    C = get_constants;
    
    
    for n = 1:length(point_dirs)
        if ~isempty(strfind(point_dirs(n).name,'.mat'))
            fn = [in_dir '/' point_dirs(n).name];
            ofn = [out_dir '/' point_dirs(n).name];
            
            if ~exist(ofn,'file') || force_transform;
            
                disp(['loading: ' fn]); tic
                surface_points = [];
                load(fn);

                if ~isempty(surface_points)

    
                    surface_points = apply_transform(T,surface_points);
%                     surface_points = apply_transform_mass(T,surface_points);
    
                    save(ofn,'surface_points');
                    toc

                end
            end
        end
       
    end
end
