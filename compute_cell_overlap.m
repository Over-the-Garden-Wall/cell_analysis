function ol_mat = compute_cell_overlap(cell_ids)
    
    num_cells = length(cell_ids);
    ol_mat = zeros(num_cells);
    
    cells = cell(num_cells,1);
    
    for n = 1:num_cells
        cells{n} = cell_data(cell_ids(n));
    end
    
    for m = 1:num_cells
        for n = m+1:num_cells
            
            sparse_mat1 = cells{m}.get_point_count_grid;
            sparse_mat2 = cells{n}.get_point_count_grid;
            
            new_size = min([size(sparse_mat1); size(sparse_mat2)]);
            
            [i,j,s] = find(sparse_mat1);
            is_good = i <= new_size(1) & j <= new_size(2);            
            sparse_mat1 = sparse(i(is_good),j(is_good),s(is_good),new_size(1),new_size(2));
            
            [i,j,s] = find(sparse_mat2);
            is_good = i <= new_size(1) & j <= new_size(2);            
            sparse_mat2 = sparse(i(is_good),j(is_good),s(is_good),new_size(1),new_size(2));
            
            ol_mat(m,n) = sum(sum(sparse_mat1.*sparse_mat2));
        end
    end
    
    if num_cells == 2
        ol_mat = ol_mat(1,2);
    else
        ol_mat = ol_mat + ol_mat';
    end
end
        