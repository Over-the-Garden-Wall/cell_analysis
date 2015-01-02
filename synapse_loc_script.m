
C = get_constants;

load('conns_synapses.mat');
conns = double(conns);

all_cell_nums = conns(1:2,:);

all_cell_nums = unique(all_cell_nums(:));

num_cells = length(all_cell_nums);

num2loc = sparse(all_cell_nums, ones(num_cells,1),1:num_cells);
soma_loc = zeros(num_cells, 3);

for c = 1:num_cells
    c_d = cell_data(all_cell_nums(c));
    
    soma_loc(c,:) = c_d.get_midpoint(true);    
end


num_synapses = size(conns,2);
thetas = zeros(num_synapses,1);
distances = zeros(num_synapses,2);

for c = 1:size(conns,2)
    rel_loc1 = (soma_loc(num2loc(conns(1,c)),2:3) - conns(5:6,c)');
    rel_loc2 = (soma_loc(num2loc(conns(2,c)),2:3) - conns(5:6,c)');
    
    distances(c,1) = sqrt(sum( rel_loc1.^2 ));
    distances(c,2) = sqrt(sum( rel_loc2.^2 ));
    
    thetas(c) = sum(rel_loc1 .* rel_loc2) / distances(c,1) / distances(c,2);
    thetas(c) = acos(thetas(c));
end
    