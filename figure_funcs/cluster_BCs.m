function [BC_partitian] = cluster_BCs

    C = get_constants;
    BCs = 60001:60230;
    
    for n = length(BCs):-1:1
        if ~exist([C.point_dir 'cell_' num2str(BCs(n)) '_surface.mat'], 'file')
            BCs(n) = [];
        end
    end
    
    num_BCs = length(BCs);
    
    %need 10th, 25th, 75, percentiles and axonal volume
    split_vals = zeros(num_BCs,4);
    
    for n = 1:num_BCs
        cell_dat = cell_data(BCs(n));
        p = cell_dat.get_surface;
        d = C.f(p(:,1));
        
        d(d<0) = [];
        
        d = sort(d);
        num_points = length(d);
        split_vals(n,1:3) = d(ceil([.1 .25 .75]*num_points));
        if cell_dat.V == 0;
            cell_dat = cell_data(BCs(n), true);
            if cell_dat.V == 0
                num_points;
            end
        end
        
        split_vals(n,4) = cell_dat.V;
    end
    
    cell_group = zeros(num_BCs,1);
    cell_group(split_vals(:,3)<36) = 12;
    cell_group(split_vals(:,3)>=36) = 34;
    
    
    is_12 = cell_group == 12;
    hist(split_vals(is_12,3)-split_vals(is_12,2), 1:20);
    
    is_1 = (split_vals(:,3)-split_vals(:,2) < 11) & is_12;
    is_2 = (split_vals(:,3)-split_vals(:,2) >= 11) & is_12;
    
    cell_group(is_1) = 1;
    cell_group(is_2) = 2;
    
    is_34 = cell_group == 34;
    
    is_4 = (split_vals(:,1) < 24) & is_34;
    is_3 = (split_vals(:,1) >= 24) & is_34;
    
    
    is_3a = (split_vals(:,4) >= 9.09*10^5) & is_3;
    is_3b = (split_vals(:,4) < 9.09*10^5) & is_3;
    
    
    
    
    
    
    
    cell_group(is_3a) = 31;
    cell_group(is_3b) = 32;
    cell_group(is_4) = 4;
    
    BC_partitian{1} = BCs(cell_group==1);
    BC_partitian{2} = BCs(cell_group==2);
    BC_partitian{3} = BCs(cell_group==31);
    BC_partitian{4} = BCs(cell_group==32);
    BC_partitian{5} = BCs(cell_group==4);
    
end
    
        
        