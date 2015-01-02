function est_length = estimate_skeleton_length(vol_voxels, surf_voxels, res)


    est_SAperVox = mean(res)^2;
    
    delta = .1;
    
    
    vol = vol_voxels*prod(res);
    
    old_rad = 0;
    est_rad = Inf;
            
    
    
    while abs(est_rad-old_rad) > delta
    
        SA = (surf_voxels*est_SAperVox);

        old_rad = est_rad;
        est_rad = 2*vol/SA;
        est_length = SA^2/4/pi/vol;
                
        
        yz_SApv = get_length_per_pix(est_rad, res([2 3])) * res(1);
        xy_SApv = get_length_per_pix(est_rad, res([1 2])) * res(3);
        xz_SApv = get_length_per_pix(est_rad, res([1 3])) * res(2);
        
        
        est_SAperVox = est_SAperVox/2 + (yz_SApv+xy_SApv+xz_SApv)/6;
%         disp(est_SAperVox);
        
    end
    
end


function len_pp = get_length_per_pix(rad, res)
    eta = .05;    

    num_vox = 2*ceil(rad./res+1)+1;
    mesh_lims = num_vox/2.*res;

    [x y] = meshgrid(-mesh_lims(2):eta:mesh_lims(2), -mesh_lims(1):eta:mesh_lims(1));
    r = sqrt(x.^2 + y.^2);
    
    is_in = r<rad;
    is_in_vox = false(num_vox);

    for vx = 1:num_vox(1)
        mesh_x = (vx-1)/eta*res(1)+1:vx/eta*res(1);
        for vy = 1:num_vox(2)
            mesh_y = (vy-1)/eta*res(2)+1:vy/eta*res(2);
            mesh_area = length(mesh_x)*length(mesh_y);
            if sum(sum(is_in(mesh_x,mesh_y)))>mesh_area/2
                is_in_vox(vx,vy) = true;
            end
        end
    end
    
    is_surf_vox = is_in_vox(2:end-1, 2:end-1) & ... 
        (is_in_vox(2:end-1, 2:end-1) ~= is_in_vox(1:end-2, 2:end-1) | ...
        is_in_vox(2:end-1, 2:end-1) ~= is_in_vox(3:end, 2:end-1) | ...
        is_in_vox(2:end-1, 2:end-1) ~= is_in_vox(2:end-1, 1:end-2) | ...
        is_in_vox(2:end-1, 2:end-1) ~= is_in_vox(2:end-1, 3:end));

    true_len = 2*pi*rad;
    len_pp = true_len / sum(is_surf_vox(:));

end
