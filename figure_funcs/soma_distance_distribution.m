function dist = soma_distance_distribution(cell_nums, bin_size)

   C = get_constants;
   num_cells = length(cell_nums);
   
   max_dist = sqrt(sum((C.max_xy - C.min_xy).^2));
   num_bins = ceil(max_dist/bin_size);
   
   
   d_count = zeros(num_bins,num_cells);
   is_valid = false(num_bins,num_cells);
   
   
   for c = 1:num_cells
       c_d1 = cell_data(cell_nums(c));
       soma1 = c_d1.get_midpoint(true);    
       soma1 = soma1(2:3);
       max_d1 = min(abs([soma1-C.max_xy, soma1-C.min_xy]));
   
       max_bin = floor(max_d1/bin_size);
       is_valid(1:max_bin, c) = true;
       
       
       for d = c+1:num_cells
           c_d2 = cell_data(cell_nums(d));
           soma2 = c_d2.get_midpoint(true);
           soma2 = soma2(2:3);
           
           d_bin = ceil(sqrt(sum((soma2-soma1).^2))/bin_size);
       
           if is_valid(d_bin, c)
               d_count(d_bin,c) = d_count(d_bin,c)+1;
           end
       end
   end
   
   dist = [];
   for k = 1:num_bins
       denom = pi*bin_size^2* (2*k - 1) * sum(is_valid(k,:));
       if denom == 0
           break
       end
       dist(k) = sum(d_count(k,:)) / denom;
   end
   
end