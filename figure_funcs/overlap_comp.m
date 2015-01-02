% close all

C = get_constants;

cell_types = {C.type.t1, C.type.t2, C.type.t3a, C.type.t3b, C.type.t4};
type_leg = {'t1', 't2', 't3a', 't3b', 't4'};

% cell_types = {[C.type.t1, C.type.t2], [C.type.t3a, C.type.t3b, C.type.t4]};
% type_leg = {'t1&2', 't3&4'};

num_types = length(cell_types);

type_strat = zeros(length(C.strat_x),num_types);

for n = 1:num_types
    num_cells = length(cell_types{n});
    for k = 1:num_cells;
        cell_dat = cell_data(cell_types{n}(k));
        my_strat = cell_dat.stratification;
        len = length(my_strat);
        
        type_strat(1:len,n) = type_strat(1:len,n) + my_strat/num_cells;
    end            
end

[quartile_data bins full_x full_data] = get_stratification_by_dist(C.type.j,[-0.6396    0.7687], 20);

rows_to_del = full_x(1) - C.strat_x(1);
type_strat = type_strat(rows_to_del+1:end,:);

if length(full_x) > size(type_strat,1)
    full_data = full_data(1:size(type_strat,1),:);
else
    type_strat = type_strat(1:length(full_x),:);
end

overlap_mat = zeros(size(full_data,2), num_types);

for t = 1:size(full_data,2)
    for n = 1:num_types
        overlap_mat(t,n) = sum(full_data(:,t).*type_strat(:,n));
    end
end


