function [angle_dat prox_dat] = get_sac_connection_angles_script

    C = get_constants;

    cell_nums = C.type.off_sac;
    num_cells = length(cell_nums);
    
    
    num_bins = 12;
    angle_bin = pi/num_bins;
    
    
    distal_thresh = 80000;
    
    angle_dat = zeros(num_bins,4);
    prox_dat = zeros(200,1);
    
    cell_dat = cell(num_cells,1);
    for n = 1:num_cells
        cell_dat{n} = cell_data(cell_nums(n));        
    end
    
    for n = 1:num_cells
        cell_mid = cell_dat{n}.get_midpoint(true);
        
        for m = 1:num_cells
            other_mid = cell_dat{m}.get_midpoint(true);
            h = sqrt(sum((other_mid(2:3)-cell_mid(2:3)).^2));
            
            p = find(cell_dat{n}.contacts(1,:)==cell_nums(m));
            
            
            for k = 1:length(p)
                conts = double(cell_dat{n}.contacts(:,p(k)))';
                
                a = sqrt(sum((conts(4:5)-cell_mid(2:3)).^2));
                b = sqrt(sum((conts(4:5)-other_mid(2:3)).^2));
                
                if a <= distal_thresh && b >= distal_thresh
                    t = 2;
                elseif a <= distal_thresh
                    t = 1;
                elseif b >= distal_thresh
                    t = 4;
                else
                    t = 3;
                end
                    
                
                    theta = acos((a^2+b^2-h^2)/2/a/b);

                    bin = ceil(theta/angle_bin);
                    angle_dat(bin,t) = angle_dat(bin,t)+1;
                    
                    prox_dat(ceil(a/1000)) = prox_dat(ceil(a/1000))+1;
%                 end
                
            end
        end
    end
end