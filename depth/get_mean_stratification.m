function [norm_strat total_strat] = get_mean_stratification(cell_nums, min_depth, max_depth)
    


    

    f = @z2ipldepth;
    
    norm_strat = zeros(max_depth-min_depth+1,1);
    total_strat = zeros(max_depth-min_depth+1,1);
    
    k = 1;
    for n = 1:length(cell_nums)
        fn = ['./surface_points_trans/cell_' num2str(cell_nums(n)) '_surface.mat'];
        if exist(fn,'file')
            load(fn);
            
            depth = ceil(f(surface_points(:,1)));
            
            my_strat = zeros(max_depth-min_depth+1,1);
            
            for k = min_depth:max_depth
                my_strat(k-min_depth+1) = sum(depth==k);
            end
            
            total_strat = total_strat+my_strat;
            norm_strat = norm_strat*(k-1)/k + my_strat/sum(my_strat)/k;
            
            k = k+1;
        end
    end        
end