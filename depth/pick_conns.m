function selected_R = pick_conns(R, target_cell, connecting_cells, cell_output)
    if ~exist('cell_output','var') || isempty(cell_output)
        cell_output = true;
    end



    latter_is_target = R(2,:) == target_cell;
    former_is_target = R(1,:) == target_cell;
    
    R(1:2,latter_is_target) = R([2 1], latter_is_target);
    R = R(:, latter_is_target | former_is_target);
    
    if cell_output 
        selected_R = cell(length(connecting_cells),1);
    else
        selected_R = [];
    end
    
    for n = 1:length(connecting_cells)
        if cell_output 
            selected_R{n} = R(:,R(2,:)==connecting_cells(n));
        else
            selected_R = [selected_R R(:,R(2,:)==connecting_cells(n))];
        end
    end
    
end